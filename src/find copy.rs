use std::{borrow::Borrow, collections::HashSet};

use bio::pattern_matching::myers::Myers;

use crate::reference::Ref;

pub fn split_cys4_regions<'a>(seq: &'a [u8], reference: &Ref) -> Option<(&'a [u8], &'a [u8])> {
    let max_dist = 6;
    
    let mut myers: Myers = reference.cys4.myers.clone();

    let mut matches: Vec<_> = myers
        .find_all(seq, max_dist).collect();

    matches
        .sort_by_key(|(_, _, dist)| *dist);

    fn disjoint(m1: &(usize, usize, u8), m2: &(usize, usize, u8)) -> bool {
        let (start1, end1, _) = *m1;
        let (start2, end2, _) = *m2;

        return start2 > end1 || start1 > end2;
    }

    let mut best: Vec<(usize, usize, u8)> = vec![];
    for m@(_, _, _) in matches {
        let dis = best.clone().into_iter()
            .map(|n| disjoint(&m, &n))
            .all(|b| b);
        
        if dis {
            best.push(m.to_owned());
        }
    }

    // get them in order
    best.sort_by_key(|(start, _, _)| *start);
    
    if let &[(_, end1, _), (start2, end2, _), (start3, _, _), ..] = &best[..] {
        Some((&seq[end1+1..start2], &seq[end2+1..start3]))
    } else {
        None
    }
}

fn split_scaffold_regions<'a>(seq: &'a [u8], reference: &Ref) -> Option<(&'a [u8], &'a [u8])> {
    todo!()
}

fn match_reference(spacer_seq: &[u8], extension_seq: &[u8], nicking_seq: &[u8], reference: Ref) -> Option<String> {
    todo!()
}

fn find(seq: &[u8], reference: Ref) -> Option<String> {

    // split into the first and second cys4 regions
    let (cys4_first, cys4_second) = split_cys4_regions(seq, &reference)?;

    // split into the first and second scaffolded regions
    let (before_scaffold, after_scaffold) = split_scaffold_regions(cys4_first, &reference)?;

    // match up the parts
    let reference = match_reference(before_scaffold, after_scaffold, cys4_second, reference);

    todo!()
}