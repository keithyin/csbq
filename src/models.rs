use std::{collections::HashMap, fs};

pub trait TUnigramModel {
    fn get_prob_with_ref_called_base(&self, ref_base: u8, called_base: u8) -> f32;

    fn get_prob_with_ref_base(&self, ref_base: u8) -> f32;
}

pub struct UnigramModel {
    refcalled2prob: HashMap<String, f32>,
}

impl UnigramModel {
    pub fn new(fname: &str) -> Self {
        let refcalled2prob = fs::read_to_string(fname)
            .unwrap()
            .split("\n")
            .into_iter()
            .map(|line| {
                let (left, right) = line.trim().split_once(" ").unwrap();
                (left.to_string(), right.parse::<f32>().unwrap())
            })
            .collect::<HashMap<_, _>>();

        Self { refcalled2prob }
    }

    pub fn get_prob(&self, key: &str) -> f32 {
        *self.refcalled2prob.get(key).unwrap()
    }
}

impl TUnigramModel for UnigramModel {
    fn get_prob_with_ref_called_base(&self, ref_base: u8, called_base: u8) -> f32 {
        let key = unsafe { String::from_utf8_unchecked(vec![ref_base, called_base]) };
        self.get_prob(&key)
    }

    fn get_prob_with_ref_base(&self, ref_base: u8) -> f32 {
        let key = unsafe { String::from_utf8_unchecked(vec![ref_base]) };
        self.get_prob(&key)
    }
}

/// AA 0.9  real A -> called A prob
/// AC 0.02
/// AG 0.02
/// AT 0.02
/// A- 0.02 real A -> deletion  prob
/// A+ 0.02 real A -> insertion prob
/// CA  ..
/// ..
/// A  0.25 prior
/// C  0.25 ..
/// G  0.25 ..
/// T  0.25 ..
pub struct DefaultUnigramModel {}

impl TUnigramModel for DefaultUnigramModel {
    fn get_prob_with_ref_called_base(&self, ref_base: u8, called_base: u8) -> f32 {
        return if ref_base == called_base { 0.9} else { 0.02 };
    }

    fn get_prob_with_ref_base(&self, _ref_base: u8) -> f32 {
        0.25
    }
}
