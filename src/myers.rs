use bio::{
    alignment::Alignment,
    pattern_matching::myers::{long, Myers},
};
use itertools::Itertools;

#[derive(Clone, Debug, Eq, PartialEq, Hash)]
pub enum VarMyers {
    Short(Myers<u64>),
    Long(long::Myers<u64>),
}

impl VarMyers {
    pub fn new(seq: &[u8]) -> Self {
        if seq.len() <= 64 {
            VarMyers::Short(Myers::<u64>::new(seq))
        } else {
            VarMyers::Long(long::Myers::<u64>::new(seq))
        }
    }

    pub fn find_best_end(&self, text: &[u8]) -> (usize, u8) {
        match self {
            VarMyers::Short(myers) => myers.find_best_end(text),
            VarMyers::Long(myers) => {
                let (a, b) = myers.find_best_end(text);
                (a, b as u8)
            }
        }
    }

    pub fn find_all(&self, seq: &[u8], edit_dist: u8) -> Vec<(usize, usize, usize)> {
        match self {
            VarMyers::Short(s) => s
                .clone()
                .find_all(seq, edit_dist)
                .map(|(s, e, d)| (s, e, d as usize))
                .collect_vec(),
            VarMyers::Long(l) => l.clone().find_all(seq, edit_dist as usize).collect_vec(),
        }
    }

    pub fn find_all_lazy<'a>(&'a mut self, text: &'a [u8], max_dist: u8) -> Option<Alignment> {
        match self {
            VarMyers::Short(myers) => {
                let mut matches = myers.find_all_lazy(text, max_dist);
                let (best_end, _) = matches.by_ref().min_by_key(|&(_, dist)| dist)?;
                let mut aln = Alignment::default();
                matches.alignment_at(best_end, &mut aln);

                Some(aln)
            }
            VarMyers::Long(myers) => {
                let mut matches = myers.find_all_lazy(text, max_dist as usize);

                let (best_end, _) = matches.by_ref().min_by_key(|&(_, dist)| dist)?;
                let mut aln = Alignment::default();
                matches.alignment_at(best_end, &mut aln);

                Some(aln)
            }
        }
    }

    pub fn find_all_disjoint(&self, seq: &[u8], edit_dist: u8) -> Vec<(usize, usize, usize)> {
        let mut matches = self.find_all(seq, edit_dist);

        matches.sort_by_key(|(_, _, dist)| *dist);

        fn disjoint(m1: &(usize, usize, usize), m2: &(usize, usize, usize)) -> bool {
            let (start1, end1, _) = *m1;
            let (start2, end2, _) = *m2;

            start2 >= end1 || start1 >= end2
        }

        let mut best: Vec<(usize, usize, usize)> = vec![];
        for m @ (_, _, _) in matches {
            let dis = best
                .clone()
                .into_iter()
                .map(|n| disjoint(&m, &n))
                .all(|b| b);

            if dis {
                best.push(m.to_owned());
            }
        }

        best
    }
}
