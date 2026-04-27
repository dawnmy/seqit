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
        }
    }

    pub fn next(&mut self) -> Result<Option<Record<'_>>> {
        self.record_buf.clear();

        if !self.has_lookahead {
            if !self.read_next_nonempty_line()? {
                return Ok(None);
            }
            self.is_fastq = match trim_crlf(&self.line_buf).first().copied() {
                Some(b'>') => false,
                Some(b'@') => true,
                _ => bail!("invalid FASTA/FASTQ record: expected '>' or '@' header"),
            };
        } else {
            std::mem::swap(&mut self.line_buf, &mut self.lookahead_line);
            self.has_lookahead = false;
        }

        let header_line = trim_crlf(&self.line_buf);
        if header_line.is_empty() {
            bail!("invalid FASTA/FASTQ record: empty header line");
        }
        self.record_buf.extend_from_slice(&header_line[1..]);
        let header_end = self.record_buf.len();

        loop {
            match self.read_next_nonempty_line_into_record_buf(true, self.is_fastq)? {
                ReadLineOutcome::Eof | ReadLineOutcome::NextHeader | ReadLineOutcome::FastqSep => {
                    break;
                }
                ReadLineOutcome::Appended(_) => {}
            }
        }
        let seq_end = self.record_buf.len();

        if self.is_fastq {
            let seq_len = seq_end - header_end;
            let mut qual_read_len = 0usize;
            while qual_read_len < seq_len {
                match self.read_next_nonempty_line_into_record_buf(false, false)? {
                    ReadLineOutcome::Eof => break,
                    ReadLineOutcome::Appended(len) => {
                        qual_read_len += len;
                        if qual_read_len > seq_len {
                            bail!(
                                "FASTQ record has sequence length {seq_len} but quality length {qual_read_len}"
                            );
                        }
                    }
                    ReadLineOutcome::NextHeader | ReadLineOutcome::FastqSep => {
                        bail!(
                            "FASTQ record has sequence length {seq_len} but quality length {qual_read_len}"
                        );
                    }
                }
            }
            if qual_read_len != seq_len {
                bail!(
                    "FASTQ record has sequence length {seq_len} but quality length {qual_read_len}"
                );
            }
        }

        let (id, desc) = parse_header(&self.record_buf[..header_end]);
        let seq = &self.record_buf[header_end..seq_end];
        let qual = self.is_fastq.then_some(&self.record_buf[seq_end..]);
        Ok(Some(Record {
            id,
            desc,
            seq,
            qual,
        }))
    }

    fn read_line_fill_buf(&mut self) -> io::Result<usize> {
        self.line_buf.clear();
        let mut total = 0usize;
        loop {
            let (consumed, done) = {
                let buf = self.reader.fill_buf()?;
                if buf.is_empty() {
                    return Ok(total);
                }
                match memchr(b'\n', buf) {
                    Some(pos) => {
                        let end = pos + 1;
                        self.line_buf.extend_from_slice(&buf[..end]);
                        (end, true)
                    }
                    None => {
                        self.line_buf.extend_from_slice(buf);
                        (buf.len(), false)
                    }
                }
            };
            self.reader.consume(consumed);
            total += consumed;
            if done {
                return Ok(total);
            }
        }
    }

    fn read_next_nonempty_line(&mut self) -> Result<bool> {
        loop {
            let n = self.read_line_fill_buf()?;
            if n == 0 {
                return Ok(false);
            }
            if !trim_crlf(&self.line_buf).is_empty() {
                return Ok(true);
            }
        }
    }

    fn read_next_nonempty_line_into_record_buf(
        &mut self,
        stop_on_fasta_header: bool,
        stop_on_fastq_sep: bool,
    ) -> Result<ReadLineOutcome> {
        loop {
            let fast_path = {
                let buf = self.reader.fill_buf()?;
                if buf.is_empty() {
                    return Ok(ReadLineOutcome::Eof);
                }
                memchr(b'\n', buf).map(|pos| {
                    let consumed = pos + 1;
                    let line = trim_crlf(&buf[..consumed]);
                    let first_char = line.first().copied();
                    let needs_copy = matches!(first_char, Some(b'>') if stop_on_fasta_header);
                    (consumed, line.len(), first_char, needs_copy)
                })
            };

            if let Some((consumed, line_len, first_char, needs_copy)) = fast_path {
                if line_len == 0 {
                    self.reader.consume(consumed);
                    continue;
                }

                if needs_copy {
                    let buf = self.reader.fill_buf()?;
                    self.lookahead_line.clear();
                    self.lookahead_line.extend_from_slice(&buf[..consumed]);
                    self.has_lookahead = true;
                    self.reader.consume(consumed);
                    return Ok(ReadLineOutcome::NextHeader);
                }

                if matches!(first_char, Some(b'+') if stop_on_fastq_sep) {
                    self.reader.consume(consumed);
                    return Ok(ReadLineOutcome::FastqSep);
                }

                let buf = self.reader.fill_buf()?;
                self.record_buf
                    .extend_from_slice(trim_crlf(&buf[..consumed]));
                self.reader.consume(consumed);
                return Ok(ReadLineOutcome::Appended(line_len));
            }

            let first_char = {
                let buf = self.reader.fill_buf()?;
                if buf.is_empty() {
                    return Ok(ReadLineOutcome::Eof);
                }
                buf[0]
            };
            if !matches!(first_char, b'>' if stop_on_fasta_header)
                && !matches!(first_char, b'+' if stop_on_fastq_sep)
            {
                return self.read_long_line_into_record_buf();
            }

            if self.read_line_fill_buf()? == 0 {
                return Ok(ReadLineOutcome::Eof);
            }
            let line = trim_crlf(&self.line_buf);
            if line.is_empty() {
                continue;
            }
            match line[0] {
                b'>' if stop_on_fasta_header => {
                    std::mem::swap(&mut self.line_buf, &mut self.lookahead_line);
                    self.has_lookahead = true;
                    return Ok(ReadLineOutcome::NextHeader);
                }
                b'+' if stop_on_fastq_sep => return Ok(ReadLineOutcome::FastqSep),
                _ => {
                    self.record_buf.extend_from_slice(line);
                    return Ok(ReadLineOutcome::Appended(line.len()));
                }
            }
        }
    }

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
                    let trimmed_len = if line.last().is_some_and(|b| *b == b'\r') {
                        line.len().saturating_sub(1)
                    } else {
                        line.len()
                    };
                    Some((pos + 1, trimmed_len, true, pos == 0, line.ends_with(b"\r")))
                } else {
                    let trimmed_len = if buf.last().is_some_and(|b| *b == b'\r') {
                        buf.len().saturating_sub(1)
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
