#[derive(Clone, Debug)]
pub struct Ran<T> {
    pub start: T,
    pub end: T,
}

impl<T> Ran<T> {
    pub fn new(start: T, end: T) -> Ran<T> {
        Ran { start, end }
    }

    pub fn map<R>(&self, f: impl Fn(&T) -> R) -> Ran<R> {
        Ran::new(f(&self.start), f(&self.end))
    }
}
