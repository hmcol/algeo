use std::cmp::Ordering;
use std::fmt::{self, Write};
use std::ops::{Add, AddAssign, Div, DivAssign, Index, Mul, MulAssign, Sub, SubAssign};

use std::collections::BTreeMap;
use std::iter::FromIterator;

use super::num::*;

// multidegree -----------------------------------------------------------------

/// type for indexing the indeterminates
type I = u8;

/// type of the degrees in a multidegree
type D = i8;

/// wrapper for `Vec<P>` to represent multidegree of monomial terms in a multivariate polynomial ring
#[derive(Clone, PartialEq, Eq, Debug, Hash)]
pub struct MDeg(pub BTreeMap<I, D>);

impl MDeg {
    /// returns empty multidegree: `[]`
    ///
    /// used for constant terms, i.e., terms without indeterminates
    pub fn zero() -> Self {
        MDeg(BTreeMap::new())
    }

    /// returns multidegree with `deg` repeated: `[deg ... deg]`
    fn repeated(deg: D, n: I) -> Self {
        MDeg((0..n).map(|i| (i, deg)).collect())
    }

    /// returns multidegree of all ones: `[1 ... 1]`
    ///
    /// used for the term `x_1` `x_2` \dots `x_n`, where each indeterminate occurs as a factor exactly once
    pub fn ones(n: I) -> Self {
        Self::repeated(1, n)
    }

    /// returns iterator over the individual degree components
    pub fn degs(&self) -> impl Iterator<Item = (&I, &D)> {
        self.0.iter()
    }

    /// returns mutable iterator over the individual degree components
    pub fn degs_mut(&mut self) -> impl Iterator<Item = (&I, &mut D)> {
        self.0.iter_mut()
    }

    /// returns the 'total degree' of a multidegree, i.e., the sum of the individual degree components
    pub fn total_deg(&self) -> D {
        self.0.values().sum()
    }

    pub fn max_idx(&self) -> I {
        *self.0.keys().max().unwrap_or(&0)
    }
}

impl AddAssign<&Self> for MDeg {
    fn add_assign(&mut self, other: &Self) {
        for (idx, other_deg) in other.degs() {
            if let Some(self_deg_ref) = self.0.get_mut(idx) {
                *self_deg_ref += other_deg;
            }
        }
    }
}

impl AddAssign for MDeg {
    fn add_assign(&mut self, other: Self) {
        *self += &other;
    }
}

impl Add for MDeg {
    type Output = Self;

    fn add(mut self, other: Self) -> Self::Output {
        self += other;

        // return
        self
    }
}

impl Add for &MDeg {
    type Output = MDeg;

    fn add(self, other: Self) -> Self::Output {
        let mut out = self.clone();

        out += other;

        // return
        out
    }
}

// term ------------------------------------------------------------------------

#[derive(Clone, PartialEq, Eq, Debug, Hash)]
pub struct Term<F: Field> {
    coef: F,
    mdeg: MDeg,
}

impl<F: Field> Term<F> {
    pub fn from_coef_mdeg(coef: F, mdeg: MDeg) -> Self {
        Term { coef, mdeg }
    }

    /// returns a term with the given coefficient and multidegree `0`
    pub fn constant(coef: F) -> Self {
        Term::from_coef_mdeg(coef, MDeg::zero())
    }

    /// hmmm 0 should not have a multidegree so this feels wrong
    pub fn zero() -> Self {
        Term::constant(F::ZERO)
    }

    // returns the constant term `1`
    pub fn one() -> Self {
        Term::constant(F::ONE)
    }

    /// returns the term of given multidegree, with coefficient `1`
    pub fn monic(mdeg: MDeg) -> Self {
        Term::from_coef_mdeg(F::ONE, mdeg)
    }

    /// maps the polynomial element `self` of `F[X]` to the corresponding polynomial function `F -> F`, given by evaluation, which is then evaluated at `x`
    pub fn eval(&self, x: &[F]) -> F {
        if x.len() < self.mdeg.max_idx().into() {
            panic!("incorrectly sized argument {:?} passed to term {:?}", x, self);
        }

        self.mdeg
            .degs()
            .map(|(&idx, &deg)| x[idx as usize].powi32(deg.into()))
            .fold(self.coef, |acc, x| acc * x)
    }
}

impl<F: Field> MulAssign for Term<F> {
    fn mul_assign(&mut self, other: Self) {
        *self *= &other;
    }
}

impl<F: Field> MulAssign<&Self> for Term<F> {
    fn mul_assign(&mut self, other: &Self) {
        self.coef *= other.coef;
        self.mdeg += &other.mdeg;
    }
}

impl<F: Field> Mul for Term<F> {
    type Output = Self;

    fn mul(mut self, other: Self) -> Self::Output {
        self *= other;

        // return
        self
    }
}


// polynomial ------------------------------------------------------------------

/// a polynomial over N variables with coefficients in the field F
#[derive(Clone, Debug, Hash)]
pub struct Polynomial<F: Field> {
    terms: Vec<Term<F>>,
}

impl<F: Field> Polynomial<F> {
    pub fn from_vec(v: Vec<Term<F>>) -> Self {
        Polynomial { terms: v }
    }

    /// returns a polynomial with only the term `t`
    pub fn monomial(t: Term<F>) -> Self {
        Self::from_vec(vec![t])
    }

    /// returns a polynomial with no terms
    pub fn zero() -> Self {
        Self::from_vec(Vec::new())
    }

    /// returns the polynomial with only the single term `1`, which correctly has the degree `[0 ... 0]`
    pub fn one() -> Self {
        Self::monomial(Term::one())
    }

    pub fn term_count(&self) -> usize {
        self.terms.len()
    }

    pub fn iter_terms(&self) -> impl Iterator<Item = &Term<F>> {
        self.terms.iter()
    }

    pub fn eval(&self, x: &[F]) -> F {
        self.terms
            .iter()
            .map(|t| t.eval(x))
            .fold(F::ZERO, |acc, y| acc + y) // should impl Sum for Field
    }
}

/// the following addition operations are implemented:
///
/// poly += term
/// poly += poly
/// poly += &poly
/// poly + poly
/// &poly + &poly

impl<F: Field> AddAssign<Term<F>> for Polynomial<F> {
    fn add_assign(&mut self, rhs: Term<F>) {
        *self += &rhs;
    }
}

impl<F: Field> AddAssign<&Term<F>> for Polynomial<F> {
    fn add_assign(&mut self, rhs: &Term<F>) {
        if let Some(term_ref) = self.terms.iter_mut().find(|t| t.mdeg == rhs.mdeg) {
            term_ref.coef += rhs.coef; // consider adding checks for zero coef
        } else {
            self.terms.push(rhs.clone());
        }
    }
}

impl<F: Field> AddAssign for Polynomial<F> {
    fn add_assign(&mut self, other: Self) {
        for t in other.terms {
            *self += t
        }
    }
}

impl<F: Field> AddAssign<&Self> for Polynomial<F> {
    fn add_assign(&mut self, other: &Self) {
        for t in other.iter_terms() {
            *self += t;
        }
    }
}

impl<F: Field> Add for Polynomial<F> {
    type Output = Self;

    fn add(mut self, other: Self) -> Self::Output {
        self += other;

        // return
        self
    }
}

impl<F: Field> Add for &Polynomial<F> {
    type Output = Polynomial<F>;

    fn add(self, other: Self) -> Self::Output {
        let mut out = self.clone();
        out += other;

        // return
        out
    }
}

/// the following multiplication operations are implemented:
///
/// poly *= term
/// poly * term
/// &poly * term
/// &poly * &poly
/// poly * poly

impl<F: Field> MulAssign<Term<F>> for Polynomial<F> {
    fn mul_assign(&mut self, rhs: Term<F>) {
        *self *= &rhs;
    }
}

impl<F: Field> MulAssign<&Term<F>> for Polynomial<F> {
    fn mul_assign(&mut self, rhs: &Term<F>) {
        for term_ref in self.terms.iter_mut() {
            *term_ref *= rhs;
        }
    }
}

impl<F: Field> Mul<&Term<F>> for Polynomial<F> {
    type Output = Self;

    fn mul(mut self, rhs: &Term<F>) -> Self::Output {
        self *= rhs;

        // return
        self
    }
}

impl<F: Field> Mul<Term<F>> for Polynomial<F> {
    type Output = Self;

    fn mul(mut self, rhs: Term<F>) -> Self::Output {
        self * &rhs
    }
}

impl<F: Field> Mul<&Term<F>> for &Polynomial<F> {
    type Output = Polynomial<F>;

    fn mul(self, rhs: &Term<F>) -> Self::Output {
        let mut out = self.clone();

        out *= rhs;

        // return
        out
    }
}

impl<F: Field> Mul<Term<F>> for &Polynomial<F> {
    type Output = Polynomial<F>;

    fn mul(self, rhs: Term<F>) -> Self::Output {
        self * &rhs
    }
}

impl<F: Field> Mul for &Polynomial<F> {
    type Output = Polynomial<F>;

    #[allow(clippy::suspicious_arithmetic_impl)]
    fn mul(self, other: Self) -> Self::Output {
        other
            .terms
            .iter()
            .fold(Polynomial::zero(), |acc, t| acc + self * t)
    }
}

impl<F: Field> Mul for Polynomial<F> {
    type Output = Self;

    fn mul(self, other: Self) -> Self::Output {
        &self * &other
    }
}