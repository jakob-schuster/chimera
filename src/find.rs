use std::collections::HashMap;

use clap::error;
use itertools::Itertools;
use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};
use crate::reference::{EfficientGuide, EfficientGuides, FinalGuides, Mismatch, Pattern, Ref};

pub fn split_cys4_regions<'a>(
    seq: &'a [u8], 
    reference: &Ref,
    error_rate: f32
) -> Option<(&'a [u8], &'a [u8])> {
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

pub fn match_reference_all_carefully<'a>(
    spacer_seq: &[u8], 
    extension_seq: &[u8], 
    nicking_seq: &[u8],
    efficient_guides: &'a EfficientGuides,
    error_rate: f32
) -> Vec<(&'a str, Mismatch)> {
    // println!("handling seqs {}, {}, {}", _seq_to_string(&spacer_seq), _seq_to_string(&extension_seq), _seq_to_string(&nicking_seq));

    fn all_matching_names<'a>(
        seq: &[u8], 
        patterns: &'a [EfficientGuide],
        error_rate: f32
    ) -> HashMap<&'a str, Mismatch> {

        let tolerance = 0.1;

        let all_matches = patterns.iter()
            .filter_map(|EfficientGuide { pattern, names }|
                pattern.clone().get_matches(seq, error_rate)
                    .into_iter()
                    .min_by(|(_, _, dist), (_, _, b_dist)| dist.partial_cmp(b_dist).unwrap())
                    .map(|a| (names, a)))
                    .sorted_by(|(_,(_, _, dist)), (_,(_, _, b_dist))| dist.partial_cmp(b_dist).unwrap());
                    // .sorted_by_key(|(_, (_, _, dist))| *dist);
        
        if let Some((_, (_, _, best_dist))) = all_matches.clone().next() {
            all_matches
                .take_while(|(_, (_, _, dist))| dist.error_rate() <= best_dist.error_rate() + tolerance)
                .flat_map(|(v, (_, _, dist))| v.iter().map(|n| (&n[..], dist)).collect_vec())
                .collect()
        } else {
            HashMap::new()
        }
    }

    let best_spacer_names = 
        all_matching_names(spacer_seq, &efficient_guides.spacers, error_rate);
    let best_extension_names = 
        all_matching_names(extension_seq, &efficient_guides.extensions, error_rate);
    let best_nicking_names = 
        all_matching_names(nicking_seq, &efficient_guides.nickings, error_rate);

    let mut final_names = HashMap::new();
    for (spacer_name, Mismatch { len: spacer_len, dist: spacer_dist }) in best_spacer_names {
        if let Some(Mismatch { len: extension_len, dist: extension_dist }) = best_extension_names.get(&spacer_name) {
            if let Some(Mismatch { len: nicking_len, dist: nicking_dist } ) = best_nicking_names.get(&spacer_name) {

                final_names.insert(spacer_name, Mismatch::new(spacer_len + extension_len + nicking_len, spacer_dist + extension_dist + nicking_dist));
            }
        }
    }

    if let Some((_, best_dist)) = final_names.clone().into_iter().min_by(|(_,dist), (_, b_dist)| dist.partial_cmp(b_dist).unwrap()) {
        final_names.into_iter().filter(|(_, dist)| *dist <= best_dist).collect_vec()
    } else {
        Vec::new()
    }
}

fn best_matches(matches: &[(usize, usize, Mismatch)], tolerance: &f32) -> Vec<(usize, usize, Mismatch)> {
    if let Some((_, _, best)) = matches.to_owned().iter()
        .min_by(|(_, _, d1), (_, _, d2)| d1.partial_cmp(d2).unwrap()) {
        
        matches.iter()
            .filter(|(_, _, dist)| dist.error_rate() <= best.error_rate() + tolerance)
            .map(|a| a.to_owned()).collect_vec()
    } else {
        vec![]
    }
}

pub fn match_reference_all_quickly<'a>(
    spacer_seq: &[u8], 
    extension_seq: &[u8], 
    nicking_seq: &[u8],
    guides: &'a FinalGuides,
    error_rate: f32
) -> Vec<(&'a str, Mismatch)> {


    fn f_a<'a>(
        name_mismatch: &(&'a str, Mismatch), 
        guides: &HashMap<String, Pattern>, 
        seq: &[u8], 
        error_rate: f32
    ) -> Option<(&'a str, Mismatch)> {
        let (name, mismatch) = name_mismatch;
        
        let pattern = guides.get(*name)?;
        
        match &best_matches(&pattern.get_matches(seq, error_rate), &0.0)[..] {
            [] => None,
            [(_, _, first), ..] => Some((
                name, 
                Mismatch { 
                    len: mismatch.len + first.len, 
                    dist: mismatch.dist + first.dist 
                }
            )),
        }
    }

    fn f_b<'a>(
        name_mismatch: &(&'a str, Mismatch), 
        guides: &HashMap<String, Pattern>, 
        seq: &[u8], 
        error_rate: f32
    ) -> Option<(&'a str, Mismatch)> {
        let (name, mismatch) = name_mismatch;
        
        let pattern = guides.get(*name)?;
        
        pattern.get_best_match(seq, error_rate).map(|m| (
            *name,
            Mismatch {
                len: mismatch.len + m.len,
                dist: mismatch.dist + m.dist,
            }
        ))
    }


    // all adequate spacers
    // let spacer_names: Vec<_> = guides.spacers.keys()
    //     .map(|name| (&name[..], Mismatch::new(0, 0)))
    //     .flat_map(|name| f_a(&name, &guides.spacers, spacer_seq, error_rate*0.5)).collect();

    // all names of sequences
    let final_names: Vec<_> = guides.spacers.keys()
        .map(|name| (&name[..], Mismatch::new(0, 0)))
        .flat_map(|name| f_b(&name, &guides.spacers, spacer_seq, error_rate))
        .flat_map(|name_mismatch| 
            f_b(&name_mismatch, &guides.extensions, extension_seq, error_rate))
        .flat_map(|name_mismatch| 
            f_b(&name_mismatch, &guides.nickings, nicking_seq, error_rate)).collect();

    if let Some((_, best)) = final_names.iter().min_by_key(|(_, dist)| dist) {
        final_names.iter().filter(|(_, dist)| dist <= best)
            .map(|(name, mismatch)| (*name, mismatch.to_owned())).collect()
    } else {
        vec![]
    }
}

#[derive(Clone, PartialEq, Eq, Hash, Debug, Default)]
pub enum StructureResult<'a> {
    // (.. scaf .. cys4 .. scaf ..) structure can be found
    WellStructured(RefResult<'a>),

    // no structure can be found
    #[default]
    BadlyStructured,
}

#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub enum RefResult<'a> {
    // there is exactly one best (spacer, extension, nicking) matching a valid triple from the reference
    Valid(&'a str, Mismatch),

    // there are no best (spacer, extension, nicking) matching a triple from the reference
    Chimera,
    
    // there are multiple best (spacer, extension, nicking) which match different reference triples
    Ambiguous,
}

pub fn structure_classify_carefully<'a>(
    seq: &[u8],
    reference: &Ref,
    efficient_guides: &'a EfficientGuides,
    error_rate: f32
) -> StructureResult<'a> {
    match break_into_regions(seq, reference, error_rate) {
        Some((spacer_seq, extension_seq, nicking_seq)) => 
            StructureResult::WellStructured(reference_classify_carefully(
                spacer_seq, extension_seq, nicking_seq, efficient_guides, error_rate
            )),
        None =>
            StructureResult::BadlyStructured,
    }
}

pub fn reference_classify_carefully<'a>(
    spacer_seq: &[u8], 
    extension_seq: &[u8], 
    nicking_seq: &[u8], 
    efficient_guides: &'a EfficientGuides,
    error_rate: f32
) -> RefResult<'a> {
    match &match_reference_all_carefully(spacer_seq, extension_seq, nicking_seq, efficient_guides, error_rate)[..] {
        [] => RefResult::Chimera,
        [(name, error)] => RefResult::Valid(name, error.clone()),
        _ => RefResult::Ambiguous
    }
}

pub fn structure_classify_quickly<'a>(
    seq: &[u8],
    reference: &Ref,
    final_guides: &'a FinalGuides,
    error_rate: f32
) -> StructureResult<'a> {
    match break_into_regions(seq, reference, error_rate) {
        Some((spacer_seq, extension_seq, nicking_seq)) => 
            StructureResult::WellStructured(reference_classify_quickly(
                spacer_seq, extension_seq, nicking_seq, final_guides, error_rate
            )),
        None =>
            StructureResult::BadlyStructured,
    }
}

pub fn reference_classify_quickly<'a>(
    spacer_seq: &[u8], 
    extension_seq: &[u8], 
    nicking_seq: &[u8], 
    final_guides: &'a FinalGuides,
    error_rate: f32
) -> RefResult<'a> {
    match &match_reference_all_quickly(spacer_seq, extension_seq, nicking_seq, final_guides, error_rate)[..] {
        [] => RefResult::Chimera,
        [(name, error)] => RefResult::Valid(*name, *error),
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
    let (extension_seq, nicking_seq) = split_cys4_regions(after_scaffold, reference, error_rate)?;
    
    // take the before-scaffold part, and split that too
    let (_, spacer_seq) = split_cys4_regions(before_scaffold, reference, error_rate)?;

    Some((spacer_seq, extension_seq, nicking_seq))
}

pub fn break_into_regions_relaxed<'a>(seq: &'a [u8], reference: &Ref, error_rate: f32) -> Option<(&'a [u8], &'a [u8], &'a [u8])> {
    // first, find scaffolds. there should be two
    let (spacer_seq, after_scaffold) = split_scaffold_regions(seq, reference, error_rate)?;

    // take the after-scaffold part, and split it by the location of cys4
    let (extension_seq, nicking_seq) = split_cys4_regions(after_scaffold, reference, error_rate)?;    

    Some((spacer_seq, extension_seq, nicking_seq))
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