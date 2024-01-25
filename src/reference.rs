use std::{collections::{HashSet, HashMap}, hash::Hash, borrow::Borrow};

use bio::pattern_matching::myers::{Myers, long};
use itertools::Itertools;

use crate::Args;

#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct Mismatch {
    pub len: usize,
    pub dist: usize,
}

impl Mismatch {
    pub fn new(len: usize, dist: usize) -> Self {
        Mismatch { len, dist }
    }

    pub fn error_rate(&self) -> f32 {
        (self.dist as f32) / (self.len as f32)
    }
}

#[derive(Clone, PartialEq, Eq, Hash)]
pub enum VarMyers {
    Short(Myers::<u64>),
    Long(long::Myers::<u64>)
}

impl VarMyers {
    fn new(seq: &[u8]) -> Self {
        if seq.len() <= 64 {
            VarMyers::Short(Myers::<u64>::new(seq))
        } else {
            VarMyers::Long(long::Myers::<u64>::new(seq))
        }
    }
    

    fn find_all(&mut self, seq: &[u8], edit_dist: u8) -> Vec<(usize, usize, u8)> {        
        match self {
            VarMyers::Short(s) => s.find_all(seq, edit_dist).collect_vec(),
            VarMyers::Long(l) => l.find_all(seq, edit_dist as usize)
                .map(|(s, e, d)| (s, e, d as u8)).collect_vec(),
        }
    }
}

#[derive(Clone, PartialEq, Eq, Hash)]
pub struct Pattern {
    pub seq: Vec<u8>,
    pub myers: VarMyers
}

#[derive(Clone)]
pub struct NamedPattern {
    pub name: String,
    pub pattern: Pattern
}

pub struct Ref {
    pub cys4: Pattern,
    pub scaffold: Pattern,
    pub guides: Vec<Guide>
}

#[derive(Clone)]
pub struct Guide {
    pub name: String,
    pub spacer: Pattern,
    pub extension: Pattern,
    pub nicking: Pattern
}

impl Guide {
    fn new(record: TSVRecord) -> Self {
        Guide {
            name: record.name,
            spacer: Pattern::new(record.spacer.as_bytes(), ),
            extension: Pattern::new(record.extension.as_bytes()),
            nicking: Pattern::new(record.nicking.as_bytes()),
        }
    }
}

#[derive(Clone)]
pub struct EfficientGuide {
    pub pattern: Pattern,
    pub names: Vec<String>
}

#[derive(Clone)]
pub struct EfficientGuides {
    pub spacers: Vec<EfficientGuide>,
    pub extensions: Vec<EfficientGuide>,
    pub nickings: Vec<EfficientGuide>,
}

impl NamedPattern {
    fn new(name: &str, seq: &[u8]) -> NamedPattern {
        NamedPattern { 
            name: String::from(name),
            pattern: Pattern::new(seq) 
        }
    }
}

impl Pattern {
    fn new(seq: &[u8]) -> Pattern {
        Pattern {
            seq: seq.to_owned(),
            myers: VarMyers::new(seq),
        }
    }

    pub fn get_matches(&self, seq: &[u8], error_rate: f32) -> Vec<(usize, usize, Mismatch)> {
        let edit_dist = (error_rate * (self.seq.len() as f32)).floor() as u8;

        let mut matches = self.myers.clone().find_all(seq, edit_dist).into_iter()
            .map(|(s, e, d)| (s, e, Mismatch::new(self.seq.len(), d as usize)))
            .collect_vec();
        matches
            .sort_by_key(|(_, _, dist)| dist.clone());

        fn disjoint<T>(m1: &(usize, usize, T), m2: &(usize, usize, T)) -> bool {
            let (start1, end1, _) = *m1;
            let (start2, end2, _) = *m2;

            start2 > end1 || start1 > end2
        }

        let mut best: Vec<(usize, usize, Mismatch)> = vec![];
        for m@(_, _, _) in matches {
            let dis = best.clone().into_iter()
                .map(|n| disjoint(&m, &n))
                .all(|b| b);
            
            if dis {
                best.push(m.clone());
            }
        }

        best
    }
}

impl Ref {
    pub fn new(arg: &Args) -> Ref {
        Ref {
            cys4: Pattern::new(arg.cys4.as_bytes()),
            scaffold: Pattern::new(arg.scaffold.as_bytes()),
            guides: parse_reference(&arg.reference_tsv),
        }
    }
}

#[derive(Debug, serde::Deserialize)]
struct TSVRecord {
    name: String,
    spacer: String,
    extension: String,
    nicking: String,
}

fn parse_reference(reference_tsv: &str) -> Vec<Guide> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_path(reference_tsv)
        .expect("Bad reference file!");

    let mut guides = Vec::new();

    for result in reader.deserialize() {
        // read a record
        let record: TSVRecord = result
            .expect("Bad reference row!");

        // add it to the guides
        guides.push(Guide::new(record));
    }

    guides
}

impl EfficientGuides {
    pub fn new(guides: &[Guide]) -> EfficientGuides {
        fn make_guides(guides: &[Guide], f: fn(&Guide) -> Pattern) -> Vec<EfficientGuide> {
            let all_names = guides.iter().map(|g| (g.name.clone(), f(g)));

            all_names.group_by(|(_, pattern)| pattern.to_owned()).into_iter()
                .map(|(pattern, group)| 
                    EfficientGuide { 
                        pattern: pattern.clone(), 
                        names: group.map(|(name, _)| name).collect_vec()})
                .collect_vec()
        }

        EfficientGuides { 
            spacers: make_guides(guides, |g| g.spacer.clone()), 
            extensions: make_guides(guides, |g| g.extension.clone()), 
            nickings: make_guides(guides, |g| g.nicking.clone())
        }
    }
}

pub struct FinalGuides {
    pub spacers: HashMap<String, Pattern>,
    pub extensions: HashMap<String, Pattern>,
    pub nickings: HashMap<String, Pattern>,
}

impl FinalGuides {
    pub fn new(guides: &[Guide]) -> FinalGuides {
        FinalGuides {
            spacers: guides.iter()
            .map(|Guide {name, spacer, extension: _, nicking: _ }| 
                (name.clone(), spacer.clone())
            ).collect(),
            extensions: guides.iter()
            .map(|Guide {name, spacer: _, extension, nicking: _ }| 
                (name.clone(), extension.clone())
            ).collect(),
            nickings: guides.iter()
            .map(|Guide {name, spacer: _, extension: _, nicking }| 
                (name.clone(), nicking.clone())
            ).collect(),
        }
    }
}