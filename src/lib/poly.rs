use std::ops::{Add, Sub, Mul, Div};
use std::collections::BTreeMap;

use super::num::*;


type P = i32;


#[derive(Clone, Debug)]
pub struct MultiDegree<const D: usize>([P; D]);

impl<const D: usize> MultiDegree<D> {
    fn total_degree(&self) -> i32 {
        0
    }
}


impl<const D: usize> Add for MultiDegree<D> {
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        let mut pows = self.0.clone();

        for (k, &b) in other.0.iter() {
            if let Some(a) = pows.get_mut(k) {
                *a += b;
            } else {
                pows.insert(*k, b);
            }
        }

        MultiDegree::new(pows)
    }
}

impl Add for &MultiDegree {
    type Output = MultiDegree;

    fn add(self, other: Self) -> Self::Output {
        let mut pows = self.pows.clone();

        for (k, &b) in other.pows.iter() {
            if let Some(a) = pows.get_mut(k) {
                *a += b;
            } else {
                pows.insert(*k, b);
            }
        }

        MultiDegree::new(pows)
    }
}


#[derive(Clone, Debug)]
pub struct Term<F: Field> {
    coef: F,
    mdeg: MultiDegree,
}

impl<F: Field> Term<F> {
    pub fn from_coef_mdeg(coef: F, mdeg: MultiDegree) -> Self {
        Term {
            coef,
            mdeg,
        }
    }

}

impl<F: Field> Mul for Term<F> {
    type Output = Self;

    fn mul(self, other: Self) -> Self::Output {
        Term::from_coef_mdeg(
            self.coef * other.coef,
            self.mdeg + other.mdeg,
        )
    }
}