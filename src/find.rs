use std::collections::HashMap;

use clap::error;
use itertools::Itertools;
use crate::reference::{Ref, EfficientGuides, EfficientGuide};

pub fn split_cys4_regions<'a>(
    seq: &'a [u8], 
    reference: &Ref,
    error_rate: f32
) -> Option<(&'a [u8], &'a [u8])> {
    let max_dist = 6;
    
    let mut matches = reference.cys4.clone().get_matches(seq, error_rate);

    // get them in order
    matches.sort_by_key(|(start, _, _)| *start);
    
    match &matches[..] {
        &[(start, end, _)] => 
            Some((&seq[..start], &seq[end..])),
        _ => None
    }
}

fn split_scaffold_regions<'a>(
    seq: &'a [u8], reference: &Ref, error_rate: f32
) -> Option<(&'a [u8], &'a [u8])> {
    let mut matches = reference.scaffold.clone().get_matches(seq, error_rate);

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

// pub fn match_reference_all(
//     spacer_seq: &[u8], 
//     extension_seq: &[u8], 
//     nicking_seq: &[u8], 
//     efficient_guides: &EfficientGuides,
//     error_rate: f32
// ) -> Vec<String> {
//     // println!("handling seqs {}, {}, {}", _seq_to_string(&spacer_seq), _seq_to_string(&extension_seq), _seq_to_string(&nicking_seq));

//     fn all_matching_names(
//         seq: &[u8], 
//         patterns: Vec<EfficientGuide>,
//         error_rate: f32
//     ) -> Vec<String> {

//         let tolerance = 0;

//         let all_matches = patterns.into_iter()
//             .map(|EfficientGuide { pattern, names }|
//                 pattern.clone().get_matches(seq, error_rate)
//                     .into_iter()
//                     .min_by_key(|(_, _, dist)| *dist)
//                     .and_then(|a| Some((names, a))))
//                     .filter_map(|o| o)
//                     .sorted_by_key(|(_, (_, _, dist))| *dist);
        
//         if let Some((_, (_, _, best_dist))) = all_matches.clone().next() {
//             all_matches
//                 .take_while(|(_, (_, _, dist))| *dist <= best_dist + tolerance)
//                 .map(|(v, _)| v)
//                 .concat()
//         } else {
//             Vec::new()
//         }
//     }

//     let best_spacer_names = 
//         all_matching_names(spacer_seq, efficient_guides.spacers.clone(), error_rate);
//     let best_extension_names = 
//         all_matching_names(extension_seq, efficient_guides.extensions.clone(), error_rate);
//     let best_nicking_names = 
//         all_matching_names(nicking_seq, efficient_guides.nickings.clone(), error_rate);

//     let names = best_spacer_names.into_iter()
//         .filter(|name| best_extension_names.contains(name))
//         .filter(|name| best_nicking_names.contains(name))
//         .collect_vec();
    
//     names
// }

pub fn match_reference_all_carefully(
    spacer_seq: &[u8], 
    extension_seq: &[u8], 
    nicking_seq: &[u8],
    efficient_guides: &EfficientGuides,
    error_rate: f32
) -> Vec<(String, f32)> {
    // println!("handling seqs {}, {}, {}", _seq_to_string(&spacer_seq), _seq_to_string(&extension_seq), _seq_to_string(&nicking_seq));

    fn all_matching_names(
        seq: &[u8], 
        patterns: Vec<EfficientGuide>,
        error_rate: f32
    ) -> HashMap<String, f32> {

        let tolerance = 0.1;

        let all_matches = patterns.into_iter()
            .map(|EfficientGuide { pattern, names }|
                pattern.clone().get_matches(seq, error_rate)
                    .into_iter()
                    .min_by(|(_, _, dist), (_, _, b_dist)| dist.partial_cmp(b_dist).unwrap())
                    .and_then(|a| Some((names, a))))
                    .filter_map(|o| o)
                    .sorted_by(|(_,(_, _, dist)), (_,(_, _, b_dist))| dist.partial_cmp(b_dist).unwrap());
                    // .sorted_by_key(|(_, (_, _, dist))| *dist);
        
        if let Some((_, (_, _, best_dist))) = all_matches.clone().next() {
            all_matches
                .take_while(|(_, (_, _, dist))| *dist <= best_dist + tolerance)
                .map(|(v, (_, _, dist))| v.into_iter().map(|n| (n, dist)).collect_vec())
                .flatten()
                .collect()
        } else {
            HashMap::new()
        }
    }

    let best_spacer_names = 
        all_matching_names(spacer_seq, efficient_guides.spacers.clone(), error_rate);
    let best_extension_names = 
        all_matching_names(extension_seq, efficient_guides.extensions.clone(), error_rate);
    let best_nicking_names = 
        all_matching_names(nicking_seq, efficient_guides.nickings.clone(), error_rate);

    let mut final_names = HashMap::new();
    for (spacer_name, spacer_dist) in best_spacer_names {
        if let Some(extension_dist) = best_extension_names.get(&spacer_name) {
            if let Some(nicking_dist) = best_nicking_names.get(&spacer_name) {
                final_names.insert(spacer_name, spacer_dist + extension_dist + nicking_dist);
            }
        }
    }

    if let Some((_, best_dist)) = final_names.clone().into_iter().min_by(|(_,dist), (_, b_dist)| dist.partial_cmp(b_dist).unwrap()) {
        final_names.into_iter().filter(|(_, dist)| dist <= &best_dist).collect_vec()
    } else {
        Vec::new()
    }
}

#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub enum StructureResult {
    // (.. scaf .. cys4 .. scaf ..) structure can be found
    WellStructured(RefResult),

    // no structure can be found
    BadlyStructured,
}

#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub enum RefResult {
    // there is exactly one best (spacer, extension, nicking) matching a valid triple from the reference
    Valid(String,u8),

    // there are no best (spacer, extension, nicking) matching a triple from the reference
    Chimera,
    
    // there are multiple best (spacer, extension, nicking) which match different reference triples
    Ambiguous,
}

pub fn structure_classify(
    seq: &[u8],
    reference: &Ref,
    efficient_guides: &EfficientGuides,
    error_rate: f32
) -> StructureResult {
    match break_into_regions(seq, reference, error_rate) {
        Some((spacer_seq, extension_seq, nicking_seq)) => 
            StructureResult::WellStructured(reference_classify(
                spacer_seq, extension_seq, nicking_seq, efficient_guides, error_rate
            )),
        None =>
            StructureResult::BadlyStructured,
    }
}

pub fn reference_classify(
    spacer_seq: &[u8], 
    extension_seq: &[u8], 
    nicking_seq: &[u8], 
    efficient_guides: &EfficientGuides,
    error_rate: f32
) -> RefResult {
    match &match_reference_all_carefully(spacer_seq, extension_seq, nicking_seq, efficient_guides, error_rate)[..] {
        [] => RefResult::Chimera,
        [(name, error)] => RefResult::Valid(name.clone(), (*error * 100.0) as u8),
        _ => RefResult::Ambiguous
    }
}

/// Tries as best as possible to detect chimeras
pub fn chimeric(
    spacer_seq: &[u8], 
    extension_seq: &[u8], 
    nicking_seq: &[u8], 
    efficient_guides: &EfficientGuides,
    error_rate: f32
) -> bool {
    // a read is chimeric if it isn't a good match with any of the references
    match_reference_all_carefully(spacer_seq, extension_seq, nicking_seq, efficient_guides, error_rate).is_empty()
}

pub fn break_into_regions<'a>(seq: &'a [u8], reference: &Ref, error_rate: f32) -> Option<(&'a [u8], &'a [u8], &'a [u8])> {
    // first, find scaffolds. there should be two
    let (before_scaffold, after_scaffold) = split_scaffold_regions(seq, reference, error_rate)?;

    // take the after-scaffold part, and split it by the location of cys4
    let (cys4_first, cys4_second) = split_cys4_regions(after_scaffold, reference, error_rate)?;

    Some((before_scaffold, cys4_first, cys4_second))
}

/// Uses a shortcut. All positives are true; some negatives are false. ie, there was a superior alignment 
pub fn quick_chimeric(
    seq: &[u8],
    reference: &Ref,
    error_rate: f32
) -> Option<bool> {
    // first, find scaffolds. there should be two
    let (before_scaffold, after_scaffold) = split_scaffold_regions(seq, reference, error_rate)?;

    // take the after-scaffold part, and split it by the location of cys4
    let (cys4_first, cys4_second) = split_cys4_regions(after_scaffold, reference, error_rate)?;

    todo!()
}