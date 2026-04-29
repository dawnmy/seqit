// Fast buffered FASTA/FASTQ reader adapted from the parser strategy in
// https://github.com/shenwei356/fastx. It keeps one reusable record buffer and
// returns borrowed record slices so streaming commands avoid per-record copies.
use anyhow::{bail, Result};
use memchr::{memchr, memchr2};
use std::io::{self, BufRead};

pub(super) struct Reader<R: BufRead> {
    reader: R,
    is_fastq: bool,
    record_buf: Vec<u8>,
    line_buf: Vec<u8>,
    lookahead_line: Vec<u8>,
    has_lookahead: bool,
    parse_id: bool,
}

pub(super) struct Record<'a> {
    pub id: &'a [u8],
    pub desc: Option<&'a [u8]>,
    pub seq: &'a [u8],
    pub qual: Option<&'a [u8]>,
}

impl<R: BufRead> Reader<R> {
    pub fn from_reader(reader: R) -> Self {
        Self {
            reader,
            is_fastq: false,
            record_buf: Vec::with_capacity(1 << 20),
            line_buf: Vec::with_capacity(1024),
            lookahead_line: Vec::with_capacity(1024),
            has_lookahead: false,
            parse_id: true,
        }
    }

    pub fn skip_id_parsing(&mut self) {
        self.parse_id = false;
    }

    pub fn next(&mut self) -> Result<Option<Record<'_>>> {
        self.record_buf.clear();

        if !self.has_lookahead {
            if !self.read_next_nonempty_line()? {
                return Ok(None);
            }
            self.is_fastq = match self.line_buf.first().copied() {
                Some(b'>') => false,
                Some(b'@') => true,
                _ => bail!("invalid FASTA/FASTQ record: expected '>' or '@' header"),
            };
        } else {
            std::mem::swap(&mut self.line_buf, &mut self.lookahead_line);
            self.has_lookahead = false;
        }

        if self.line_buf.is_empty() {
            bail!("invalid FASTA/FASTQ record: empty header line");
        }
        self.record_buf.extend_from_slice(&self.line_buf[1..]);
        let header_end = self.record_buf.len();

        if self.is_fastq {
            loop {
                match self.read_next_nonempty_line_into_record_buf::<true, true>()? {
                    ReadLineOutcome::Eof
                    | ReadLineOutcome::NextHeader
                    | ReadLineOutcome::FastqSep => break,
                    ReadLineOutcome::Appended(_) => {}
                }
            }
        } else {
            loop {
                match self.read_next_nonempty_line_into_record_buf::<true, false>()? {
                    ReadLineOutcome::Eof
                    | ReadLineOutcome::NextHeader
                    | ReadLineOutcome::FastqSep => break,
                    ReadLineOutcome::Appended(_) => {}
                }
            }
        }
        let seq_end = self.record_buf.len();

        if self.is_fastq {
            let seq_len = seq_end - header_end;
            let mut qual_read_len = 0usize;
            while qual_read_len < seq_len {
                match self.read_qual_line_into_record_buf()? {
                    ReadLineOutcome::Eof => break,
                    ReadLineOutcome::Appended(len) => {
                        qual_read_len += len;
                        if qual_read_len > seq_len {
                            bail!(
                                "FASTQ record has sequence length {seq_len} but quality length {qual_read_len}"
                            );
                        }
                    }
                    ReadLineOutcome::NextHeader | ReadLineOutcome::FastqSep => unreachable!(),
                }
            }
            if qual_read_len != seq_len {
                bail!(
                    "FASTQ record has sequence length {seq_len} but quality length {qual_read_len}"
                );
            }
        }

        let (id, desc) = if self.parse_id {
            parse_header(&self.record_buf[..header_end])
        } else {
            (&self.record_buf[..header_end], None)
        };
        let seq = &self.record_buf[header_end..seq_end];
        let qual = self.is_fastq.then_some(&self.record_buf[seq_end..]);
        Ok(Some(Record {
            id,
            desc,
            seq,
            qual,
        }))
    }

    #[inline(always)]
    fn read_line_fill_buf(&mut self) -> io::Result<usize> {
        self.line_buf.clear();

        let mut total = 0usize;
        loop {
            let buf = self.reader.fill_buf()?;
            if buf.is_empty() {
                if self.line_buf.last() == Some(&b'\r') {
                    self.line_buf.pop();
                }
                return Ok(total);
            }

            let (consumed, done) = match memchr(b'\n', buf) {
                Some(pos) => {
                    let end = pos + 1;
                    self.line_buf.extend_from_slice(&buf[..end]);
                    (end, true)
                }
                None => {
                    self.line_buf.extend_from_slice(buf);
                    (buf.len(), false)
                }
            };

            self.reader.consume(consumed);
            total += consumed;

            if done {
                self.line_buf.pop();
                if self.line_buf.last() == Some(&b'\r') {
                    self.line_buf.pop();
                }
                return Ok(total);
            }
        }
    }

    #[inline(always)]
    fn read_next_nonempty_line(&mut self) -> Result<bool> {
        loop {
            let n = self.read_line_fill_buf()?;
            if n == 0 {
                return Ok(false);
            }
            if !self.line_buf.is_empty() {
                return Ok(true);
            }
        }
    }

    #[inline(always)]
    fn read_next_nonempty_line_into_record_buf<
        const STOP_ON_FASTA_HEADER: bool,
        const STOP_ON_FASTQ_SEP: bool,
    >(
        &mut self,
    ) -> Result<ReadLineOutcome> {
        loop {
            let buf = self.reader.fill_buf()?;
            if buf.is_empty() {
                return Ok(ReadLineOutcome::Eof);
            }

            if let Some(pos) = memchr(b'\n', buf) {
                let consumed = pos + 1;
                let line_len = trim_crlf(&buf[..consumed]).len();

                if line_len == 0 {
                    self.reader.consume(consumed);
                    continue;
                }

                let first_char = buf[0];
                if STOP_ON_FASTA_HEADER && first_char == b'>' {
                    self.lookahead_line.clear();
                    self.lookahead_line.extend_from_slice(&buf[..line_len]);
                    self.reader.consume(consumed);
                    self.has_lookahead = true;
                    return Ok(ReadLineOutcome::NextHeader);
                }

                if STOP_ON_FASTQ_SEP && first_char == b'+' {
                    self.reader.consume(consumed);
                    return Ok(ReadLineOutcome::FastqSep);
                }

                self.record_buf.extend_from_slice(&buf[..line_len]);
                self.reader.consume(consumed);
                return Ok(ReadLineOutcome::Appended(line_len));
            }

            let first_char = buf[0];
            if !((STOP_ON_FASTA_HEADER && first_char == b'>')
                || (STOP_ON_FASTQ_SEP && first_char == b'+'))
            {
                return self.read_long_line_into_record_buf();
            }

            if self.read_line_fill_buf()? == 0 {
                return Ok(ReadLineOutcome::Eof);
            }
            if self.line_buf.is_empty() {
                continue;
            }

            let first_char = self.line_buf[0];
            if STOP_ON_FASTA_HEADER && first_char == b'>' {
                std::mem::swap(&mut self.line_buf, &mut self.lookahead_line);
                self.has_lookahead = true;
                return Ok(ReadLineOutcome::NextHeader);
            }
            if STOP_ON_FASTQ_SEP && first_char == b'+' {
                return Ok(ReadLineOutcome::FastqSep);
            }

            let len = self.line_buf.len();
            self.record_buf.extend_from_slice(&self.line_buf);
            return Ok(ReadLineOutcome::Appended(len));
        }
    }

    #[inline(always)]
    fn read_qual_line_into_record_buf(&mut self) -> Result<ReadLineOutcome> {
        loop {
            let buf = self.reader.fill_buf()?;
            if buf.is_empty() {
                return Ok(ReadLineOutcome::Eof);
            }

            if let Some(pos) = memchr(b'\n', buf) {
                let consumed = pos + 1;
                let line_len = trim_crlf(&buf[..consumed]).len();
                if line_len == 0 {
                    self.reader.consume(consumed);
                    continue;
                }
                self.record_buf.extend_from_slice(&buf[..line_len]);
                self.reader.consume(consumed);
                return Ok(ReadLineOutcome::Appended(line_len));
            }

            return self.read_long_line_into_record_buf();
        }
    }

    #[inline(always)]
    fn read_long_line_into_record_buf(&mut self) -> Result<ReadLineOutcome> {
        let mut line_len = 0usize;
        let mut pending_cr = false;

        loop {
            let step = {
                let buf = self.reader.fill_buf()?;
                if buf.is_empty() {
                    if line_len == 0 && !pending_cr {
                        return Ok(ReadLineOutcome::Eof);
                    }
                    None
                } else if let Some(pos) = memchr(b'\n', buf) {
                    let line = &buf[..pos];
                    let trimmed_len = if line.ends_with(b"\r") {
                        line.len() - 1
                    } else {
                        line.len()
                    };
                    Some((pos + 1, trimmed_len, true, pos == 0, line.ends_with(b"\r")))
                } else {
                    let trimmed_len = if buf.ends_with(b"\r") {
                        buf.len() - 1
                    } else {
                        buf.len()
                    };
                    Some((buf.len(), trimmed_len, false, false, buf.ends_with(b"\r")))
                }
            };

            let Some((consumed, trimmed_len, has_lf, lf_only, ends_with_cr)) = step else {
                return Ok(ReadLineOutcome::Appended(line_len));
            };

            let buf = self.reader.fill_buf()?;
            let data = if has_lf {
                &buf[..consumed - 1]
            } else {
                &buf[..consumed]
            };
            let data = &data[..trimmed_len];

            if pending_cr && !(has_lf && lf_only) {
                self.record_buf.push(b'\r');
                line_len += 1;
            }
            if !data.is_empty() {
                self.record_buf.extend_from_slice(data);
                line_len += data.len();
            }

            pending_cr = ends_with_cr && !has_lf;
            self.reader.consume(consumed);

            if has_lf {
                if line_len == 0 {
                    continue;
                }
                return Ok(ReadLineOutcome::Appended(line_len));
            }
        }
    }
}

enum ReadLineOutcome {
    Eof,
    Appended(usize),
    NextHeader,
    FastqSep,
}

#[inline]
fn trim_crlf(line: &[u8]) -> &[u8] {
    let mut end = line.len();
    if end > 0 && line[end - 1] == b'\n' {
        end -= 1;
    }
    if end > 0 && line[end - 1] == b'\r' {
        end -= 1;
    }
    &line[..end]
}

#[inline]
fn parse_header(line: &[u8]) -> (&[u8], Option<&[u8]>) {
    let Some(id_end) = memchr2(b' ', b'\t', line) else {
        return (line, None);
    };
    let mut desc_start = id_end;
    while desc_start < line.len() && matches!(line[desc_start], b' ' | b'\t') {
        desc_start += 1;
    }
    let desc = (desc_start < line.len()).then_some(&line[desc_start..]);
    (&line[..id_end], desc)
}

#[cfg(test)]
mod tests {
    use super::Reader;
    use std::io::{BufReader, Cursor};

    type OwnedRecord = (Vec<u8>, Option<Vec<u8>>, Vec<u8>, Option<Vec<u8>>);

    fn read_to_owned(input: &[u8], capacity: usize) -> anyhow::Result<Vec<OwnedRecord>> {
        let cursor = Cursor::new(input);
        let mut reader = Reader::from_reader(BufReader::with_capacity(capacity, cursor));
        let mut out = Vec::new();
        while let Some(record) = reader.next()? {
            out.push((
                record.id.to_vec(),
                record.desc.map(<[u8]>::to_vec),
                record.seq.to_vec(),
                record.qual.map(<[u8]>::to_vec),
            ));
        }
        Ok(out)
    }

    #[test]
    fn reads_wrapped_fastq_with_small_buffer_and_no_final_lf() {
        let records = read_to_owned(b"@r1 desc\nAC\nGT\n+\nII\nII", 3).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].0, b"r1");
        assert_eq!(records[0].1.as_deref(), Some(&b"desc"[..]));
        assert_eq!(records[0].2, b"ACGT");
        assert_eq!(records[0].3.as_deref(), Some(&b"IIII"[..]));
    }

    #[test]
    fn can_skip_id_parsing_for_raw_header_output() {
        let cursor = Cursor::new(b">r1   desc\tmore\nACGT\n");
        let mut reader = Reader::from_reader(BufReader::with_capacity(8, cursor));
        reader.skip_id_parsing();

        let record = reader.next().unwrap().unwrap();
        assert_eq!(record.id, b"r1   desc\tmore");
        assert_eq!(record.desc, None);
        assert_eq!(record.seq, b"ACGT");
    }

    #[test]
    fn keeps_line_buffer_small_for_long_sequence_lines() {
        let seq = b"ACGT".repeat(4096);
        let mut input = b">r1\n".to_vec();
        input.extend_from_slice(&seq);
        input.extend_from_slice(b"\n>r2\nTGCA");
        let cursor = Cursor::new(input);
        let mut reader = Reader::from_reader(BufReader::with_capacity(3, cursor));

        let first = reader.next().unwrap().unwrap();
        assert_eq!(first.id, b"r1");
        assert_eq!(first.seq, seq);
        assert!(reader.line_buf.capacity() < seq.len() / 2);

        let second = reader.next().unwrap().unwrap();
        assert_eq!(second.id, b"r2");
        assert_eq!(second.seq, b"TGCA");
    }

    #[test]
    fn accepts_crlf_and_skips_blank_lines() {
        let records = read_to_owned(b"\r\n>r1\r\nAC\r\n\r\nGT\r\n>r2\r\nTGCA\r\n", 4).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].2, b"ACGT");
        assert_eq!(records[1].2, b"TGCA");
    }

    #[test]
    fn rejects_nonempty_invalid_first_line() {
        let err = read_to_owned(b" \n", 8).unwrap_err().to_string();
        assert!(err.contains("invalid FASTA/FASTQ record"));
    }

    #[test]
    fn rejects_unequal_fastq_lengths() {
        let err = read_to_owned(b"@r1\nACGT\n+\nIII\n", 8)
            .unwrap_err()
            .to_string();
        assert!(err.contains("sequence length 4 but quality length 3"));
    }
}
