use itertools::Itertools;
use crate::reference::{Ref, Guide, Pattern};

pub fn split_cys4_regions<'a>(seq: &'a [u8], reference: &Ref) -> Option<(&'a [u8], &'a [u8])> {
    let max_dist = 6;
    
    let mut matches = reference.cys4.clone().get_matches(seq, max_dist);

    // get them in order
    matches.sort_by_key(|(start, _, _)| *start);
    
    match &matches[..] {
        &[(start, end, _)] => 
            Some((&seq[..start], &seq[end..])),
        _ => None
    }
}

fn split_scaffold_regions<'a>(
    seq: &'a [u8], reference: &Ref
) -> Option<(&'a [u8], &'a [u8])> {
    let max_dist = 20;
    
    let mut matches = reference.scaffold.clone().get_matches(seq, max_dist);

    // get them in order
    matches.sort_by_key(|(start, _, _)| *start);
    
    match &matches[..] {
        &[(start1, end1, _), (start2, _, _)] => 
            Some((&seq[..start1], &seq[end1+1..start2])),
        // also catch cases where there was just one
        // &[(start1, end1, _)] => 
        //     Some((&seq[..start1], &seq[end1+1..])),
        _ => None
    }
}

/// Tries as best as possible to match the reference
fn match_reference(
    spacer_seq: &[u8], 
    extension_seq: &[u8], 
    nicking_seq: &[u8], 
    reference: &Ref
) -> Option<String> {
    let mut name = None;
    for guide in reference.guides.clone() {
        if !guide.spacer.clone().get_matches(spacer_seq, 4).is_empty() &&
            !guide.extension.clone().get_matches(extension_seq, 4).is_empty() &&
            !guide.nicking.clone().get_matches(nicking_seq, 4).is_empty() {
            
            // return one the minute you find it
            name = Some(String::from(guide.name));
            break
        }
    }
    
    name
}

/// Tries as best as possible to detect chimeras
pub fn chimeric(
    spacer_seq: &[u8], 
    extension_seq: &[u8], 
    nicking_seq: &[u8], 
    reference: &Ref
) -> bool {
    fn contains_match(seq: &[u8], patterns: Vec<(String, Pattern)>) -> Vec<String> {
        let tolerance = 5;

        let a = patterns.into_iter()
            .map(| (name, pattern) | pattern.clone()
                .get_matches(seq, 10).into_iter()
                    .min_by_key(|(_, _, dist)| *dist).and_then(|a| Some((name, a))))
            .filter_map(|o| o)
            .sorted_by_key(|(_, (_, _, dist))| *dist);

        if let Some((_, (_, _, best_dist))) = a.clone().next() {
            a
                .filter(|(_, (_, _, dist))| *dist <= best_dist + tolerance)
                .map(|(n, _)| n)
                .collect_vec()
        } else {
            Vec::new()
        }
    }

    let best_spacer_names = contains_match(spacer_seq, reference.guides.clone()
        .into_iter()
            .map(| Guide { name, spacer, extension, nicking } | 
                (name, spacer)).collect_vec());

    let best_extension_names = contains_match(extension_seq, reference.guides.clone()
        .into_iter()
            .map(| Guide { name, spacer, extension, nicking } | 
                (name, extension)).collect_vec());

    let best_nicking_names = contains_match(spacer_seq, reference.guides.clone()
        .into_iter()
            .map(| Guide { name, spacer, extension, nicking } | 
                (name, nicking)).collect_vec());

    let names = best_spacer_names.into_iter()
        .filter(|name| best_extension_names.contains(name))
        .filter(|name| best_nicking_names.contains(name))
        .collect_vec();
    
    names.is_empty()

    // let mut found_match = false;
    // for spacer_name in &best_spacer_names {
    //     for extension_name in &best_extension_names {
    //         for nicking_name in &best_nicking_names {
    //             if [spacer_name, extension_name, nicking_name].into_iter().all_equal() {
    //                 found_match = true;
    //             }
    //         }
    //     }
    // }

    // !found_match
    // let results = [best_spacer_names, best_extension_names, best_nicking_names];
    
}

pub fn break_into_regions<'a>(seq: &'a [u8], reference: &Ref) -> Option<(&'a [u8], &'a [u8], &'a [u8])> {
    // first, find scaffolds. there should be two
    let (before_scaffold, after_scaffold) = split_scaffold_regions(seq, reference)?;

    // take the after-scaffold part, and split it by the location of cys4
    let (cys4_first, cys4_second) = split_cys4_regions(after_scaffold, reference)?;

    Some((before_scaffold, cys4_first, cys4_second))
}

pub fn find(seq: &[u8], reference: &Ref) -> Option<String> {
    // first divide it up into regions
    let (before_scaffold, cys4_first, cys4_second) = break_into_regions(seq, reference)?;

    // try to match against all the references
    match_reference(before_scaffold, cys4_first, cys4_second, reference)
}