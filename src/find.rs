use itertools::Itertools;
use crate::reference::{Ref, Guide, Pattern, EfficientGuides, EfficientGuide};

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

pub fn match_reference_all(
    spacer_seq: &[u8], 
    extension_seq: &[u8], 
    nicking_seq: &[u8], 
    efficient_guides: &EfficientGuides
) -> Vec<String> {
    fn all_matching_names(
        seq: &[u8], 
        patterns: Vec<EfficientGuide>
    ) -> Vec<String> {
        let tolerance = 0;

        let all_matches = patterns.into_iter()
            .map(|EfficientGuide { pattern, names }|
                pattern.clone().get_matches(seq, 10)
                    .into_iter()
                    .min_by_key(|(_, _, dist)| *dist)
                    .and_then(|a| Some((names, a))))
                    .filter_map(|o| o)
                    .sorted_by_key(|(_, (_, _, dist))| *dist);
        
        if let Some((_, (_, _, best_dist))) = all_matches.clone().next() {
            all_matches
                .take_while(|(_, (_, _, dist))| *dist <= best_dist + tolerance)
                .map(|(v, _)| v)
                .concat()
        } else {
            Vec::new()
        }
    }

    let best_spacer_names = 
        all_matching_names(spacer_seq, efficient_guides.spacers.clone());
    let best_extension_names = 
        all_matching_names(extension_seq, efficient_guides.extensions.clone());
    let best_nicking_names = 
        all_matching_names(nicking_seq, efficient_guides.nickings.clone());

    let names = best_spacer_names.into_iter()
        .filter(|name| best_extension_names.contains(name))
        .filter(|name| best_nicking_names.contains(name))
        .collect_vec();
    
    names
}

#[derive(Clone)]
pub enum StructureResult {
    // (scaf, cys4, scaf) structure can be found
    WellStructured(RefResult),

    // no structure can be found
    BadlyStructured,
}

#[derive(Clone)]
pub enum RefResult {
    // there is exactly one best (spacer, extension, nicking) matching a valid triple from the reference
    Valid(String),

    // there are no best (spacer, extension, nicking) matching a triple from the reference
    Chimera,
    
    // there are multiple best (spacer, extension, nicking) which match different reference triples
    Ambiguous,
}

pub fn structure_classify(
    seq: &[u8],
    reference: &Ref,
    efficient_guides: &EfficientGuides
) -> StructureResult {
    match break_into_regions(seq, reference) {
        Some((spacer_seq, extension_seq, nicking_seq)) => 
            StructureResult::WellStructured(reference_classify(
                spacer_seq, extension_seq, nicking_seq, efficient_guides
            )),
        None =>
            StructureResult::BadlyStructured,
    }
}

pub fn reference_classify(
    spacer_seq: &[u8], 
    extension_seq: &[u8], 
    nicking_seq: &[u8], 
    efficient_guides: &EfficientGuides
) -> RefResult {
    match &match_reference_all(spacer_seq, extension_seq, nicking_seq, efficient_guides)[..] {
        [] => RefResult::Chimera,
        [name] => RefResult::Valid(name.clone()),
        _ => RefResult::Ambiguous
    }
}

/// Tries as best as possible to detect chimeras
pub fn chimeric(
    spacer_seq: &[u8], 
    extension_seq: &[u8], 
    nicking_seq: &[u8], 
    efficient_guides: &EfficientGuides
) -> bool {
    // a read is chimeric if it isn't a good match with any of the references
    match_reference_all(spacer_seq, extension_seq, nicking_seq, efficient_guides).is_empty()
}

pub fn break_into_regions<'a>(seq: &'a [u8], reference: &Ref) -> Option<(&'a [u8], &'a [u8], &'a [u8])> {
    // first, find scaffolds. there should be two
    let (before_scaffold, after_scaffold) = split_scaffold_regions(seq, reference)?;

    // take the after-scaffold part, and split it by the location of cys4
    let (cys4_first, cys4_second) = split_cys4_regions(after_scaffold, reference)?;

    Some((before_scaffold, cys4_first, cys4_second))
}
