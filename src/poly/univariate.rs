use std::ops::{Add, Mul, Neg, Sub};

use crate::core::num::{One, Zero};

trait Ring:
    Sized
    + Copy
    + PartialEq
    + Add<Output = Self>
    + Sub<Output = Self>
    + Neg<Output = Self>
    + Zero
    + Mul<Output = Self>
    + One
{
}

#[derive(Debug, Clone, Copy, PartialEq, Hash)]
struct Term<R: Ring> {
    coef: R,
    deg: u8,
}

impl<R: Ring> Term<R> {
    fn new(coef: R, deg: u8) -> Self {
        Self { coef, deg }
    }

    fn zero_of_deg(deg: u8) -> Self {
        Self::new(R::ZERO, deg)
    }
}

struct UnivariatePolynomial<R: Ring> {
    terms: Vec<R>,
}

impl<R: Ring> UnivariatePolynomial<R> {
    fn new(terms: Vec<R>) -> Self {
        Self { terms }
    }
}

impl<R: Ring> Zero for UnivariatePolynomial<R> {
    const ZERO: Self = Self { terms: vec![] };

    fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }
}

impl<R: Ring> Add for UnivariatePolynomial<R> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        let len = self.terms.len().max(rhs.terms.len());
        let mut terms = Vec::with_capacity(len);

        for d in 0..len {
            let lhs_term = self.terms.get(d).copied().unwrap_or(R::ZERO);
            let rhs_term = rhs.terms.get(d).copied().unwrap_or(R::ZERO);

            terms.push(lhs_term + rhs_term);
        }

        if let Some(deg) = terms.iter().rposition(R::is_zero) {
            terms.truncate(deg);
            Self::new(terms)
        } else {
            Self::ZERO
        }
    }
}

impl<R: Ring> Neg for UnivariatePolynomial<R> {
    type Output = Self;

    fn neg(mut self) -> Self {
        self.terms.iter_mut().for_each(|x| *x = -*x);
        self
    }
}

impl<R: Ring> Sub for UnivariatePolynomial<R> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        self + -rhs
    }
}