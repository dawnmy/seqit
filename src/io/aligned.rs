use std::alloc::{alloc_zeroed, dealloc, handle_alloc_error, Layout};
use std::io::{self, BufRead, Read};
use std::ptr::NonNull;
use std::slice;

pub(crate) struct AlignedBufReader<R> {
    inner: R,
    buf: AlignedBuffer,
    pos: usize,
    filled: usize,
}

impl<R: Read> AlignedBufReader<R> {
    pub(crate) fn with_capacity_alignment(capacity: usize, alignment: usize, inner: R) -> Self {
        Self {
            inner,
            buf: AlignedBuffer::new(capacity.max(1), alignment.max(1)),
            pos: 0,
            filled: 0,
        }
    }
}

impl<R: Read> Read for AlignedBufReader<R> {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        if out.is_empty() {
            return Ok(0);
        }
        if self.pos == self.filled && out.len() >= self.buf.len() {
            return self.inner.read(out);
        }

        let available = self.fill_buf()?;
        let n = available.len().min(out.len());
        out[..n].copy_from_slice(&available[..n]);
        self.consume(n);
        Ok(n)
    }
}

impl<R: Read> BufRead for AlignedBufReader<R> {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        if self.pos >= self.filled {
            self.filled = self.inner.read(self.buf.as_mut_slice())?;
            self.pos = 0;
        }
        Ok(&self.buf.as_slice()[self.pos..self.filled])
    }

    fn consume(&mut self, amt: usize) {
        self.pos = self.filled.min(self.pos.saturating_add(amt));
    }
}

struct AlignedBuffer {
    ptr: NonNull<u8>,
    len: usize,
    layout: Layout,
}

// SAFETY: AlignedBuffer owns its allocation and only exposes it through &self/&mut self.
unsafe impl Send for AlignedBuffer {}

impl AlignedBuffer {
    fn new(len: usize, alignment: usize) -> Self {
        let alignment = alignment.next_power_of_two();
        let layout =
            Layout::from_size_align(len, alignment).expect("valid aligned IO buffer layout");
        let ptr = unsafe {
            let ptr = alloc_zeroed(layout);
            NonNull::new(ptr).unwrap_or_else(|| handle_alloc_error(layout))
        };
        Self { ptr, len, layout }
    }

    fn len(&self) -> usize {
        self.len
    }

    fn as_slice(&self) -> &[u8] {
        unsafe { slice::from_raw_parts(self.ptr.as_ptr(), self.len) }
    }

    fn as_mut_slice(&mut self) -> &mut [u8] {
        unsafe { slice::from_raw_parts_mut(self.ptr.as_ptr(), self.len) }
    }
}

impl Drop for AlignedBuffer {
    fn drop(&mut self) {
        unsafe {
            dealloc(self.ptr.as_ptr(), self.layout);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::AlignedBufReader;
    use std::io::{BufRead, Cursor, Read};

    #[test]
    fn buffer_is_aligned_and_reads_lines() {
        let mut reader = AlignedBufReader::with_capacity_alignment(8, 64, Cursor::new(b"a\nbc"));
        let ptr = reader.fill_buf().unwrap().as_ptr() as usize;
        assert_eq!(ptr % 64, 0);

        let mut line = String::new();
        reader.read_line(&mut line).unwrap();
        assert_eq!(line, "a\n");
        line.clear();
        reader.read_line(&mut line).unwrap();
        assert_eq!(line, "bc");
    }

    #[test]
    fn large_reads_bypass_internal_buffer() {
        let mut reader = AlignedBufReader::with_capacity_alignment(4, 64, Cursor::new(b"abcdef"));
        let mut out = [0u8; 6];
        assert_eq!(reader.read(&mut out).unwrap(), 6);
        assert_eq!(&out, b"abcdef");
    }
}
