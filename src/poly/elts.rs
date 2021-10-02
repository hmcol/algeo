use std::ops::{Add, Div, Mul, Neg, Sub};
use xops::binop;

use super::mdeg::MultiDegree;
use crate::core::num::Field;

// structs ---------------------------------------------------------------------

/// single term of a multivariate polynomial with coefficients in `F`.
#[derive(Clone, Debug)]
pub struct Term<F: Field> {
    pub coef: F,
    pub mdeg: MultiDegree,
}

/// a multivariate polynomial with coefficients in the field `F`
#[derive(Clone, Debug)]
pub struct Polynomial<F: Field> {
    // pub? yeah pub.
    pub terms: Vec<Term<F>>,
}

// implementations -------------------------------------------------------------

impl<F: Field> Term<F> {
    /// Returns the term with the given coefficient and multidegree.
    ///
    /// WARNING: Does not sanitize zeros.
    #[inline]
    pub fn new_unchecked(coef: F, mdeg: MultiDegree) -> Self {
        Term { coef, mdeg }
    }

    /// Returns the term with the given coefficient and multidegree.
    ///
    /// NOTE: Sanitizes zeros.
    ///
    /// This is the most basic way one should create a new term struct.
    pub fn new(coef: F, mdeg: MultiDegree) -> Self {
        if coef == F::ZERO {
            Term::zero()
        } else {
            Term::new_unchecked(coef, mdeg)
        }
    }

    /// returns a term with the given coefficient and multidegree `MDeg::0`
    ///
    /// WARNING: Does not sanitize zeros.
    ///
    /// this is used for interpreting elements of the field as possible terms
    /// in polynomials over the field
    #[inline]
    pub fn constant_unchecked(coef: F) -> Self {
        Term::new_unchecked(coef, MultiDegree::zero())
    }

    /// Returns the constant zero term.
    #[inline]
    pub fn zero() -> Self {
        Term::constant_unchecked(F::ZERO)
    }

    /// Checks whether the term is zero.
    #[inline]
    pub fn is_zero(&self) -> bool {
        self.coef == F::ZERO
    }

    /// Returns the constant one term
    #[inline]
    pub fn one() -> Self {
        Term::constant_unchecked(F::ONE)
    }

    /// returns the monic term of given multidegree: `Term { coef: 1, mdeg }`
    ///
    /// WARNING: Does not sanitize zeros.
    pub fn monic(mdeg: MultiDegree) -> Self {
        Term::new_unchecked(F::ONE, mdeg)
    }

    /// returns a term representing a single variable/indeterminate:
    /// - `var(0, k) = x_0^k = x^k`
    /// - `var(1, k) = x_1^k = y^k`
    /// - `var(2, k) = x_2^k = z^k`
    /// - `var(j, k) = x_j^k`
    pub fn var(idx: usize, deg: u8) -> Self {
        if deg == 0 {
            Term::one()
        } else {
            Term::monic(MultiDegree::var(idx, deg))
        }
    }

    /// maps the polynomial element `self =: f ∈ F[x_1, ..., x_n]` to the
    /// corresponding
    /// polynomial function `eval_f: F^n -> F`, and the image of `x ∈ F^n`
    /// under this function is returned
    ///
    /// should be changed to return some form of `Result<F, EvaluationError>`
    pub fn eval(&self, x: &[F]) -> F {
        if x.len() < self.mdeg.len() {
            panic!(
                "incorrectly sized argument {:?} passed to term {:?}",
                x, self
            );
        }

        self.mdeg
            .degs()
            .zip(x)
            .map(|(&deg, val)| val.powi32(deg.into()))
            .fold(self.coef, Mul::mul)
    }

    pub fn divides(&self, other: &Term<F>) -> bool {
        !self.is_zero() && self.mdeg.is_succ(&other.mdeg)
    }

    pub fn try_div(&self, other: &Term<F>) -> Option<Term<F>> {
        if self.is_zero() {
            return Some(Term::zero());
        }
        if other.is_zero() {
            return None;
        }
        // `self` and `other` guaranteed nonzero

        Some(Term::new_unchecked(
            self.coef / other.coef,
            self.mdeg.checked_sub(&other.mdeg)?,
        ))
    }
}

impl<F: Field> Polynomial<F> {
    /// Returns the polynomial with the given terms.
    ///
    /// WARNING: Does not sanitize zeros.
    #[inline]
    pub fn new_unchecked(terms: Vec<Term<F>>) -> Self {
        Polynomial { terms }
    }

    /// Returns polynomial with terms from `vec`, filtered for zero coefficients
    ///
    /// NOTE: Sanitizes zeros.
    ///
    /// Ideally, this should not be necessary; polynomial computations should
    /// be careful to keep themselves clean operations should be
    /// structured such that no additional filtering is necessary except during
    /// creation
    pub fn new(terms: Vec<Term<F>>) -> Self {
        Self::new_unchecked(terms.into_iter().filter(|t| !t.is_zero()).collect())
    }

    /// Returns the zero polynomial.
    ///
    /// Currently, this is a polynomial with no terms
    #[inline]
    pub fn zero() -> Self {
        Self::new_unchecked(Vec::new())
    }

    /// should not be used until we can guarantee that zero terms get filtered
    #[inline]
    pub fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    /// Returns the polynomial with only the term `t`
    ///
    /// WARNING: Does not sanitize zeros.
    #[inline]
    pub fn term_unchecked(t: Term<F>) -> Self {
        Self::new_unchecked(vec![t])
    }

    /// Returns `coef` as a polynomial
    ///
    /// WARNING: Does not sanitize zeros.
    ///
    /// but that doesn't really mean much except for stability stuff
    #[inline]
    pub fn constant_unchecked(coef: F) -> Self {
        Self::term_unchecked(Term::constant_unchecked(coef))
    }

    pub fn var(idx: usize, deg: u8) -> Self {
        Self::term_unchecked(Term::var(idx, deg))
    }

    /// returns the polynomial with only the constant term `1`
    #[inline]
    pub fn one() -> Self {
        Self::term_unchecked(Term::one())
    }

    /// Returns an iterator over (immutably) borrowed terms: `Item = &Term<F>`.
    ///
    /// Should replace return with in-house iterator struct.
    #[inline]
    pub fn terms(&self) -> std::slice::Iter<Term<F>> {
        self.terms.iter()
    }

    /// Returns an iterator over mutably borrowed terms: `Item = &mut Term<F>`.
    ///
    /// Should replace return with in-house iterator struct.
    #[inline]
    pub fn terms_mut(&mut self) -> std::slice::IterMut<Term<F>> {
        self.terms.iter_mut()
    }

    /// Evaluates `self` in `F[x_1, ..., x_n]` at the point `x` in `F^m`.
    ///
    /// # Panics
    ///
    /// Method panics if `m < n`, i.e., if `x` does not provide values up the
    /// maximum index variable of `self`.
    #[inline]
    pub fn eval(&self, x: &[F]) -> F {
        self.terms().map(|t| t.eval(x)).sum()
    }
}

// defaults --------------------------------------------------------------------

impl_zero_default! { Term<F> where F: Field }
impl_zero_default! { Polynomial<F> where F: Field }

// conversions -----------------------------------------------------------------

impl<F: Field> From<F> for Term<F> {
    #[inline]
    fn from(coef: F) -> Self {
        Term::constant_unchecked(coef)
    }
}

impl<F: Field> From<F> for Polynomial<F> {
    fn from(coef: F) -> Self {
        if coef == F::ZERO {
            Polynomial::zero()
        } else {
            Polynomial::constant_unchecked(coef)
        }
    }
}

impl<F: Field> From<Term<F>> for Polynomial<F> {
    fn from(term: Term<F>) -> Self {
        if term.is_zero() {
            Polynomial::zero()
        } else {
            Polynomial::term_unchecked(term)
        }
    }
}

impl<F: Field> From<&Term<F>> for Polynomial<F> {
    #[inline]
    fn from(term: &Term<F>) -> Self {
        Polynomial::from(term.clone())
    }
}

// equality --------------------------------------------------------------------

impl<F: Field> PartialEq for Term<F> {
    fn eq(&self, other: &Self) -> bool {
        self.coef == other.coef && self.mdeg == other.mdeg
    }
}

impl<F: Field> Eq for Term<F> {}

// operations ------------------------------------------------------------------

// add

#[binop(derefs)]
impl<F: Field> Add for &Term<F> {
    type Output = Polynomial<F>;

    fn add(self, rhs: &Term<F>) -> Self::Output {
        if self.is_zero() {
            // could be that `rhs == 0`; must sanitize zeros
            return rhs.into();
        }
        if rhs.is_zero() {
            // could be that `self == 0`; must sanitize zeros
            return self.into();
        }
        // `self` and `rhs` guaranteed nonzero terms

        if self.mdeg == rhs.mdeg {
            let coef = self.coef + rhs.coef;

            // could be that `self.coef == -rhs.coef`; must sanitize zeros
            Term::new(coef, self.mdeg.clone()).into()
        } else {
            // know that `self` and `rhs` are nonzero with different multidegrees
            Polynomial::new_unchecked(vec![self.clone(), rhs.clone()])
        }
    }
}

#[binop(commute, derefs)]
impl<F: Field> Add<&Term<F>> for &Polynomial<F> {
    type Output = Polynomial<F>;

    fn add(self, rhs: &Term<F>) -> Self::Output {
        if self.is_zero() {
            // could be that `rhs == 0`; must sanitize zeros
            return rhs.into();
        }
        if rhs.is_zero() {
            // trust that `self` has clean zeros
            return self.clone();
        }
        // `self` and `rhs` guaranteed nonzero

        let mut terms = self.terms.clone();

        // i think this if/else is a contender for the canonical way to add a nonzero term `rhs` to a polynomial `self`
        //
        // maybe make into a method of polynomial?
        if let Some(i) = terms.iter().position(|t| t.mdeg == rhs.mdeg) {
            terms[i].coef += rhs.coef;

            // could be that `self.terms[i] == 0`; must sanitize
            if terms[i].is_zero() {
                terms.remove(i);
            }
        } else {
            // this is a clean operation:
            // - `self` has no nonzero terms of the same multidegree as `rhs`, so no simplifying terms is possible
            // - `rhs` is nonzero, so zeros remain clean
            terms.push(rhs.clone());
        }

        // trust that `self` started with clean zeros; adding `rhs` doesn't spoil
        Polynomial::new_unchecked(terms)
    }
}

#[binop(derefs)]
impl<F: Field> Add for &Polynomial<F> {
    type Output = Polynomial<F>;

    #[inline]
    fn add(self, rhs: &Polynomial<F>) -> Self::Output {
        // clean operation so long as `Poly + Term` is clean
        //
        // in fact, if both `self` and `rhs` are trusted to have clean zeros, then this can be improved to only check at the start for zeros, then fold could use a possible `add_nonzero_term` method for polynomials
        rhs.terms().fold(self.clone(), Add::add)
    }
}

// neg

impl<F: Field> Neg for Term<F> {
    type Output = Term<F>;

    fn neg(self) -> Self::Output {
        // trust that `self` has clean zeros, taking negative won't spoil
        Term::new_unchecked(-self.coef, self.mdeg)
    }
}

impl<F: Field> Neg for &Term<F> {
    type Output = Term<F>;

    fn neg(self) -> Self::Output {
        -self.clone()
    }
}

impl<F: Field> Neg for &Polynomial<F> {
    type Output = Polynomial<F>;

    #[inline]
    fn neg(self) -> Self::Output {
        // trust that `self` has clean zeros, taking negative won't spoil
        Polynomial::new_unchecked(self.terms().map(Neg::neg).collect())
    }
}

impl<F: Field> Neg for Polynomial<F> {
    type Output = Polynomial<F>;

    #[inline]
    fn neg(mut self) -> Self::Output {
        // trust that `self` has clean zeros, taking negative won't spoil
        for t in self.terms_mut() {
            t.coef = -t.coef;
        }
        self
    }
}

// sub

#[binop(derefs)]
impl<F: Field> Sub for &Term<F> {
    type Output = Polynomial<F>;

    fn sub(self, rhs: &Term<F>) -> Self::Output {
        self + &-rhs
    }
}

#[binop(derefs)]
impl<F: Field> Sub<&Polynomial<F>> for &Term<F> {
    type Output = Polynomial<F>;

    #[inline]
    fn sub(self, rhs: &Polynomial<F>) -> Self::Output {
        self + &-rhs
    }
}

#[binop(derefs)]
impl<F: Field> Sub<&Term<F>> for &Polynomial<F> {
    type Output = Polynomial<F>;

    #[inline]
    fn sub(self, rhs: &Term<F>) -> Self::Output {
        self + &-rhs
    }
}

#[binop(derefs)]
impl<F: Field> Sub for &Polynomial<F> {
    type Output = Polynomial<F>;

    #[inline]
    fn sub(self, rhs: &Polynomial<F>) -> Self::Output {
        self + &rhs.neg()
    }
}

// mul

#[binop(derefs)]
impl<F: Field> Mul for &Term<F> {
    type Output = Term<F>;

    fn mul(self, rhs: &Term<F>) -> Self::Output {
        if self.is_zero() || rhs.is_zero() {
            return Term::zero();
        }
        // `self` and `rhs` guaranteed to be nonzero
        Term::new_unchecked(self.coef * rhs.coef, &self.mdeg + &rhs.mdeg)
    }
}

#[binop(derefs)]
impl<F: Field> Mul<F> for &Polynomial<F> {
    type Output = Polynomial<F>;

    #[inline]
    fn mul(self, rhs: F) -> Self::Output {
        self * Term::constant_unchecked(rhs)
    }
}

#[binop(commute, derefs)]
impl<F: Field> Mul<&Term<F>> for &Polynomial<F> {
    type Output = Polynomial<F>;

    #[allow(clippy::suspicious_arithmetic_impl)]
    fn mul(self, rhs: &Term<F>) -> Self::Output {
        if rhs.is_zero() || self.is_zero() {
            return Polynomial::zero();
        }
        // `self` and `rhs` guaranteed nonzero

        // this is a clean operation:
        // - all terms nonzero, so zeros stay clean
        // - all multidegrees change uniformly, so no simplification is possible
        Polynomial::new_unchecked(self.terms().map(|t| t * rhs).collect())
    }
}

#[binop(derefs)]
impl<F: Field> Mul for &Polynomial<F> {
    type Output = Polynomial<F>;

    #[allow(clippy::suspicious_arithmetic_impl)]
    fn mul(self, rhs: &Polynomial<F>) -> Self::Output {
        self.terms()
            .map(|t| t * rhs)
            .reduce(Add::add)
            .unwrap_or_default()
    }
}

// div

#[binop(derefs)]
impl<F: Field> Div<F> for &Polynomial<F> {
    type Output = Polynomial<F>;

    #[allow(clippy::suspicious_arithmetic_impl)]
    #[inline]
    fn div(self, rhs: F) -> Self::Output {
        self * rhs.inv()
    }
}

// tests -----------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mdeg() {
        dbg!(MultiDegree::zero());

        /* let a = MDeg::from_pairs(&[(0, 3), (3, 5), (4, 2)]);
        let b = MDeg::from_pairs(&[(0, 1), (1, 2), (4, 1), (5, 2)]);

        let c = MDeg::from_pairs(&[(0, 4), (1, 2), (3, 5), (4, 3), (5, 2)]);

        println!("\na = {}", a);
        println!("b = {}", b);

        let d = &a + &b;

        println!("\nc = {}", d);

        assert_eq!(a + b, c); */
    }

    #[test]
    fn test_term() {
        println!("Term::zero() = {}", Term::<f64>::zero());
        println!("Term::one() = {}", Term::<f64>::one());

        let x = Term::<f64>::var(0, 1);
        let y = Term::<f64>::var(1, 1);
        let z = Term::<f64>::var(2, 1);

        println!("x = {}", x);
        println!("y = {}", y);
        println!("z = {}", z);

        let t = &(x * y) * z;

        println!("\nt = x * y * z = {}", t);

        let c = Term::from(37.0);

        println!("\nc = {}", c);

        let d = (&c * &c) * (&t * &t);

        println!("\nc^2 * t^2 = {}", d);
        println!("c^2 * t^2 = {:?}", d);
    }

    #[test]
    fn test_poly() {
        let c = |coef| Term::<f64>::constant_unchecked(coef);
        let x = |deg| Term::<f64>::var(0, deg);
        let y = |deg| Term::<f64>::var(1, deg);
        let z = |deg| Term::<f64>::var(2, deg);

        let trm = |coef, [i, j, k]: [u8; 3]| c(coef) * x(i) * y(j) * z(k);

        println!("x^4 = {}", x(4));
        println!("y^2 = {}", y(2));
        println!("z^7 = {}", z(7));

        let p = trm(5.0, [1, 2, 3]);
        let q = trm(7.0, [4, 0, 2]);
        let r = trm(13.0, [0, 5, 6]);

        println!("\np = {}", p);
        println!("q = {}", q);
        println!("r = {}", r);

        let f = Polynomial::new_unchecked(vec![p, q, r]);

        println!("\nf = {}", f);

        let p = trm(5.0, [1, 2, 3]);
        let q = trm(7.0, [4, 0, 2]);
        let r = trm(13.0, [0, 5, 6]);

        let g = p + q + r;

        println!("\nf = {}", g);
    }
}
