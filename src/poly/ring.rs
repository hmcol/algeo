use itertools::{EitherOrBoth, Itertools};
use std::cmp::Ordering;
use std::fmt;
use std::ops::{Add, Div, Mul, Neg, Sub};
use xops::binop;

use crate::core::num::*;

// structs ---------------------------------------------------------------------

/// type for indexing the indeterminates
type I = usize;

/// type of the degrees in a multidegree
///
/// may need to be changed to `u8` as many of the computations are turning out to require nonnegative degrees
type D = u8;

/// Multidegree for a monomial; wraps a `Vec<D>`.
///
/// This is the treatment of multidegrees as (nonnegative) integer tuples.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct MDeg(pub Vec<D>);

/// single term of a multivariate polynomial with coefficients in `F`.
#[derive(Clone, Debug)]
pub struct Term<F: Field> {
    pub coef: F,
    pub mdeg: MDeg,
}

/// a multivariate polynomial with coefficients in the field `F`
#[derive(Clone, Debug)]
pub struct Polynomial<F: Field> {
    // pub? yeah pub.
    pub terms: Vec<Term<F>>,
}

// implementations -------------------------------------------------------------

impl MDeg {
    /// Returns an empty multidegree.
    ///
    /// In any mathematical context, you should prefer `MDeg::zero`.
    #[inline]
    pub fn new() -> Self {
        MDeg::zero()
    }

    /// Trims trailing zeros
    #[inline]
    pub fn trim_zeros(&mut self) {
        if let Some(idx) = self.degs().rposition(|deg| *deg != 0) {
            self.0.truncate(idx + 1);
        }
    }

    /// Returns the multidegree wrapping the given vector, after truncating off
    /// any trailing zeros
    pub fn from_vec(mut vec: Vec<D>) -> Self {
        if let Some(idx) = vec.iter().rposition(|deg| *deg != 0) {
            vec.truncate(idx + 1);
            MDeg(vec)
        } else {
            MDeg::zero()
        }
    }

    /// returns the empty multidegree
    ///
    /// used for constant terms, i.e., terms without indeterminates
    #[inline]
    pub fn zero() -> Self {
        MDeg(Vec::new())
    }

    // hack until we can guarantee multidegrees never have trailing zeros
    pub fn is_zero(&self) -> bool {
        self.0.is_empty()
    }

    /// returns the multidegree with the sole entry `idx: deg`
    ///
    /// can be used for representing variables/indeterminates:
    /// - the multidegree of `x = x_0` is taken to be `{0: 1}`
    /// - the multidegree of `y = x_1` is taken to be `{1: 1}`
    /// - the multidegree of `z = x_2` is taken to be `{2: 1}`
    pub fn var(idx: I, deg: D) -> Self {
        if deg == 0 {
            MDeg::zero()
        } else {
            let mut degs = vec![0; idx];
            degs.push(deg);

            MDeg(degs)
        }
    }

    /// returns iterator over the individual degree components
    ///
    /// Should replace return with in-house iterator struct.
    #[inline]
    pub fn degs(&self) -> std::slice::Iter<D> {
        self.0.iter()
    }

    /// returns mutable iterator over the individual degree components
    ///
    /// Should replace return with in-house iterator struct.
    #[inline]
    pub fn degs_mut(&mut self) -> std::slice::IterMut<D> {
        self.0.iter_mut()
    }

    /// returns the 'total degree' of a multidegree, i.e., the sum of the
    /// individual degree components
    #[inline]
    pub fn total_deg(&self) -> D {
        self.degs().sum()
    }

    /// returns the maximum index for which `self` contains an entry.
    ///
    /// in other words, this is the minimum value `n` for which we would
    /// consider `self` to be an element of the polynomial ring in `n` variables
    #[inline]
    pub fn len(&self) -> I {
        self.0.len()
    }

    pub fn is_succ(&self, other: &MDeg) -> bool {
        if self.len() > other.len() {
            return false;
        }

        for (a, b) in self.degs().zip(other.degs()) {
            if a > b {
                return false;
            }
        }
        
        true
    }

    pub fn try_sub(&self, other: &MDeg) -> Option<MDeg> {
        if other.is_zero() {
            return Some(self.clone());
        }
        // `other` is guaranteed nonzero

        if self.is_zero() {
            // `self - other` < 0 so subtraction would fail
            return None;
        }
        // `self` and `other` are guaranteed nonzero

        let mut degs = Vec::with_capacity(self.len());
        let mut zero_cache = 0;

        for pair in self.degs().zip_longest(other.degs()) {
            let c = match pair {
                EitherOrBoth::Both(a, b) => match a.cmp(b) {
                    // `a - b < 0`; subtraction has failed
                    Ordering::Less => return None,
                    // `a - b == 0`; increment zero cache
                    Ordering::Equal => {
                        zero_cache += 1;
                        continue;
                    }
                    // `a - b > 0`; carry on with subtraction
                    Ordering::Greater => a - b,
                },
                // only `self` degs remaining, last of which is trusted to be nonzero, so just take the rest
                EitherOrBoth::Left(a) => *a,
                // only `rhs` degs remaining, last of which is trusted to be nonzero, so subtraction has failed
                EitherOrBoth::Right(_) => return None,
            };
            // `c` is guaranteed nonzero

            // catch up with zeros
            degs.append(&mut vec![0; zero_cache]);
            zero_cache = 0;

            degs.push(c);
        }

        Some(MDeg(degs))
    }
}

impl<F: Field> Term<F> {
    /// Returns the term with the given coefficient and multidegree.
    ///
    /// WARNING: Does not sanitize zeros.
    #[inline]
    pub fn new_unchecked(coef: F, mdeg: MDeg) -> Self {
        Term { coef, mdeg }
    }

    /// Returns the term with the given coefficient and multidegree.
    ///
    /// NOTE: Sanitizes zeros.
    ///
    /// This is the most basic way one should create a new term struct.
    pub fn new(coef: F, mdeg: MDeg) -> Self {
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
        Term::new_unchecked(coef, MDeg::zero())
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
    pub fn monic(mdeg: MDeg) -> Self {
        Term::new_unchecked(F::ONE, mdeg)
    }

    /// returns a term representing a single variable/indeterminate:
    /// - `var(0, k) = x_0^k = x^k`
    /// - `var(1, k) = x_1^k = y^k`
    /// - `var(2, k) = x_2^k = z^k`
    /// - `var(j, k) = x_j^k`
    pub fn var(idx: I, deg: D) -> Self {
        if deg == 0 {
            Term::one()
        } else {
            Term::monic(MDeg::var(idx, deg))
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
            self.mdeg.try_sub(&other.mdeg)?,
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
    pub fn monomial_unchecked(t: Term<F>) -> Self {
        Self::new_unchecked(vec![t])
    }

    /// Returns `coef` as a polynomial
    ///
    /// WARNING: Does not sanitize zeros.
    #[inline]
    pub fn constant_unchecked(coef: F) -> Self {
        Self::monomial_unchecked(Term::constant_unchecked(coef))
    }

    /// Returns `coef` as a polynomial
    ///
    /// NOTE: Sanitizes zeros.
    pub fn constant(coef: F) -> Self {
        if coef == F::ZERO {
            Self::zero()
        } else {
            Self::monomial_unchecked(Term::constant_unchecked(coef))
        }
    }

    /// returns the polynomial with only the constant term `1`
    #[inline]
    pub fn one() -> Self {
        Self::monomial_unchecked(Term::one())
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

/// Implements the `Default` trait by calling `Self::zero`.
///
/// For hopefully obvious reasons, this only works if the given type has an unambiguous function called `zero`.
macro_rules! impl_zero_default {
    ($Type:ty $(where $($generics:tt)*)?) => {
        impl$(<$($generics)*>)? Default for $Type {
            #[inline]
            fn default() -> Self {
                Self::zero()
            }
        }
    };
}

impl_zero_default! { MDeg }
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
            Polynomial::monomial_unchecked(term)
        }
    }
}

impl<F: Field> From<&Term<F>> for Polynomial<F> {
    fn from(term: &Term<F>) -> Self {
        if term.is_zero() {
            Polynomial::zero()
        } else {
            Polynomial::monomial_unchecked(term.clone())
        }
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

#[binop(derefs)]
impl Add for &MDeg {
    type Output = MDeg;

    fn add(self, rhs: &MDeg) -> Self::Output {
        let mut degs = Vec::with_capacity(self.len().max(rhs.len()));
        let mut zero_cache = 0;

        for pair in self.degs().zip_longest(rhs.degs()) {
            let c = match pair {
                EitherOrBoth::Both(0, 0) | EitherOrBoth::Left(0) | EitherOrBoth::Right(0) => {
                    zero_cache += 1;
                    continue;
                }
                EitherOrBoth::Both(a, b) => a + b,
                EitherOrBoth::Left(a) => *a,
                EitherOrBoth::Right(b) => *b,
            };
            // `c` is guaranteed nonzero

            // catch up with zeros
            degs.append(&mut vec![0; zero_cache]);
            zero_cache = 0;

            degs.push(c);
        }

        MDeg(degs)
    }
}

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

#[binop(refs_clone)]
impl<F: Field> Add for Term<F> {
    type Output = Polynomial<F>;

    fn add(self, rhs: Term<F>) -> Self::Output {
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
            Term::new(coef, self.mdeg).into()
        } else {
            // know that `self` and `rhs` are nonzero with different multidegrees
            Polynomial::new_unchecked(vec![self, rhs])
        }
    }
}

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

#[binop(derefs)]
impl<F: Field> Sub for &Term<F> {
    type Output = Polynomial<F>;

    fn sub(self, rhs: &Term<F>) -> Self::Output {
        self + &-rhs
    }
}

#[binop(refs_clone)]
impl<F: Field> Add<Term<F>> for Polynomial<F> {
    type Output = Polynomial<F>;

    fn add(mut self, rhs: Term<F>) -> Self::Output {
        if self.is_zero() {
            // could be that `rhs == 0`; must sanitize zeros
            return rhs.into();
        }
        if rhs.is_zero() {
            // trust that `self` has clean zeros
            return self;
        }
        // `self` and `rhs` guaranteed nonzero

        // i think this if/else is a contender for the canonical way to add a nonzero term `rhs` to a polynomial `self`
        //
        // maybe make into a method of polynomial?
        if let Some(i) = self.terms().position(|t| t.mdeg == rhs.mdeg) {
            self.terms[i].coef += rhs.coef;

            // could be that `self.terms[i] == 0`; must sanitize
            if self.terms[i].is_zero() {
                self.terms.remove(i);
            }
        } else {
            // this is a clean operation:
            // - `self` has no nonzero terms of the same multidegree as `rhs`, so no simplifying terms is possible
            // - `rhs` is nonzero, so zeros remain clean
            self.terms.push(rhs);
        }

        // trust that `self` started with clean zeros; adding `rhs` doesn't spoil
        self
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

#[binop(commute)]
impl<F: Field> Add<&Polynomial<F>> for Polynomial<F> {
    type Output = Polynomial<F>;

    #[inline]
    fn add(self, rhs: &Polynomial<F>) -> Self::Output {
        // clean operation so long as `Poly + Term` is clean
        //
        // in fact, if both `self` and `rhs` are trusted to have clean zeros, then this can be improved to only check at the start for zeros, then fold could use a possible `add_nonzero_term` method for polynomials
        rhs.terms().fold(self, Add::add)
    }
}

impl<F: Field> Add for Polynomial<F> {
    type Output = Polynomial<F>;

    #[inline]
    fn add(self, rhs: Polynomial<F>) -> Self::Output {
        self + &rhs
    }
}

impl<F: Field> Add for &Polynomial<F> {
    type Output = Polynomial<F>;

    #[inline]
    fn add(self, rhs: &Polynomial<F>) -> Self::Output {
        self.clone() + rhs
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

#[binop(commute)]
impl<F: Field> Sub<&Polynomial<F>> for Polynomial<F> {
    type Output = Polynomial<F>;

    #[inline]
    fn sub(self, rhs: &Polynomial<F>) -> Self::Output {
        self + -rhs
    }
}

impl<F: Field> Sub for Polynomial<F> {
    type Output = Polynomial<F>;

    #[inline]
    fn sub(self, rhs: Polynomial<F>) -> Self::Output {
        self - &rhs
    }
}

impl<F: Field> Sub for &Polynomial<F> {
    type Output = Polynomial<F>;

    #[inline]
    fn sub(self, rhs: &Polynomial<F>) -> Self::Output {
        self.clone() - rhs
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

#[binop(derefs)]
impl<F: Field> Mul<F> for &Polynomial<F> {
    type Output = Polynomial<F>;

    #[inline]
    fn mul(self, rhs: F) -> Self::Output {
        self * Term::constant_unchecked(rhs)
    }
}

#[binop(derefs)]
impl<F: Field> Div<F> for &Polynomial<F> {
    type Output = Polynomial<F>;

    #[allow(clippy::suspicious_arithmetic_impl)]
    #[inline]
    fn div(self, rhs: F) -> Self::Output {
        self * rhs.inv()
    }
}

// shorthands ------------------------------------------------------------------

macro_rules! fn_var_term {
    (@with_doc $var:ident, $idx:literal, $doc_str:expr) => {
        #[doc = $doc_str]
        pub fn $var<F: Field>(deg: D) -> Term<F> {
            Term::var($idx, deg)
        }
    };
    (@doc_of $var:expr, $idx:expr) => {
        concat!(
            "Shorthand for $",
            $var,
            " = x_",
            $idx,
            "\\in F[x_1, \\dots, x_n]$.\n\n",
            "As a Term, `",
            $var,
            "(d)` is monic with multidegree `{ ",
            $idx,
            ": d }`.",
            ""
        )
    };
    ($var:ident -> $idx:literal) => {
        fn_var_term! { @with_doc $var, $idx,
            fn_var_term!(@doc_of stringify!($var), stringify!($idx))
        }
    };
}

fn_var_term! { x -> 0 }
fn_var_term! { y -> 1 }
fn_var_term! { z -> 2 }

fn_var_term! { w -> 3 }
fn_var_term! { u -> 4 }
fn_var_term! { v -> 5 }

// display ---------------------------------------------------------------------

impl fmt::Display for MDeg {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.len() <= 5 {
            self.write_var(f, 0, "x")?;
            self.write_var(f, 1, "y")?;
            self.write_var(f, 2, "z")?;
            self.write_var(f, 3, "u")?;
            self.write_var(f, 4, "v")?;
            self.write_var(f, 5, "w")?;
        } else {
            write!(f, "X{:?}", self.0)?;
        }
        Ok(())
    }
}

impl MDeg {
    fn write_var(&self, f: &mut fmt::Formatter, idx: I, ident: &str) -> fmt::Result {
        if let Some(&deg) = self.0.get(idx) {
            if deg == 1 {
                write!(f, "{}", ident)?;
            } else if deg != 0 {
                write!(f, "{}{}", ident, superscript(deg))?;

                // bonus latex mode
                // write!(f, "{}^{{{}}}", ident, deg)?;
            }
        }
        Ok(())
    }
}

pub fn superscript(n: u8) -> String {
    match n {
        0 => String::from("⁰"),
        1 => String::from("¹"),
        2 => String::from("²"),
        3 => String::from("³"),
        4 => String::from("⁴"),
        5 => String::from("⁵"),
        6 => String::from("⁶"),
        7 => String::from("⁷"),
        8 => String::from("⁸"),
        9 => String::from("⁹"),
        _ => superscript(n / 10) + &superscript(n % 10),
    }
}

impl<F: Field> fmt::Display for Term<F> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.mdeg.is_zero() {
            write!(f, "{}", self.coef)
        } else if self.coef == F::ONE {
            write!(f, "{}", self.mdeg)
        } else if self.coef == -F::ONE {
            write!(f, "-{}", self.mdeg)
        } else {
            write!(f, "{}{}", self.coef, self.mdeg)
        }
    }
}

impl fmt::Display for Polynomial<f64> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut term_iter = self.terms.iter();

        if let Some(term) = term_iter.next() {
            write!(f, "{}", term)?;
        } else {
            return write!(f, "0");
        }

        for term in term_iter {
            if term.coef.is_sign_negative() {
                write!(
                    f,
                    " - {}",
                    Term::new_unchecked(-term.coef, term.mdeg.clone())
                )?;
            } else {
                write!(f, " + {}", term)?;
            }
        }

        // return
        fmt::Result::Ok(())
    }
}

impl fmt::Display for Polynomial<crate::core::frac::Frac> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut term_iter = self.terms.iter();

        if let Some(term) = term_iter.next() {
            write!(f, "{}", term)?;
        } else {
            return write!(f, "0");
        }

        for term in term_iter {
            if term.coef.numer.is_negative() {
                write!(
                    f,
                    " - {}",
                    Term::new_unchecked(-term.coef, term.mdeg.clone())
                )?;
            } else {
                write!(f, " + {}", term)?;
            }
        }

        // return
        fmt::Result::Ok(())
    }
}

// tests -----------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mdeg() {
        dbg!(MDeg::zero());

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

        let x = x::<f64>(1);
        let y = y::<f64>(1);
        let z = z::<f64>(1);

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

        let trm = |coef, [i, j, k]: [D; 3]| c(coef) * x(i) * y(j) * z(k);

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
