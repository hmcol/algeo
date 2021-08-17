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

/// wrapper for `Vec<P>` to represent multidegree of monomial terms in a
/// multivariate polynomial ring
#[derive(Clone, PartialEq, Eq, Debug, Hash)]
pub struct MDeg(pub BTreeMap<I, D>);

impl MDeg {
    /// returns the empty multidegree
    ///
    /// used for constant terms, i.e., terms without indeterminates
    pub fn zero() -> Self {
        MDeg(BTreeMap::new())
    }

    /// returns multidegree of all ones: `[1 ... 1]`
    ///
    /// used for the term `x_1` `x_2` \dots `x_n`, where each indeterminate
    /// occurs as a factor exactly once
    pub fn ones(n: I) -> Self {
        Self::repeated(1, n)
    }

    /// returns the multidegree where each of `n` entries is `deg`.
    ///
    /// explicitly, `{0: deg, ..., n-1: deg}`.
    fn repeated(deg: D, n: I) -> Self {
        MDeg((0..n).map(|i| (i, deg)).collect())
    }

    /// returns the multidegree corresponding the the given slice of tuples
    /// 
    /// if `tuples` contains the pair `(i, d)`, then the resulting multidegree
    /// will have the entry `i: d`
    ///
    /// filtering zeros ensures equality holds when the only difference between
    /// two multidegrees is the presence of some zero entries
    ///
    /// this should probably be expanded to take any iterator over tuples, and
    /// any public constructor should pass through this filter
    pub fn from_tuples(tuples: &[(I, D)]) -> Self {
        MDeg(
            tuples
                .iter()
                .copied()
                .filter(|&(_, deg)| deg != 0i8)
                .collect()
        )
    }

    /// returns the multidegree with the sole entry `idx: deg`
    /// 
    /// can be used for representing variables/indeterminates:
    /// - the multidegree of `x = x_0` is taken to be `{0: 1}`
    /// - the multidegree of `y = x_1` is taken to be `{1: 1}`
    /// - the multidegree of `z = x_2` is taken to be `{2: 1}`
    pub fn single(idx: I, deg: D) -> Self {
        MDeg::from_tuples(&[(idx, deg)])
    }

    /// returns iterator over the individual degree components
    pub fn degs(&self) -> impl Iterator<Item = (&I, &D)> {
        self.0.iter()
    }

    /// returns mutable iterator over the individual degree components
    pub fn degs_mut(&mut self) -> impl Iterator<Item = (&I, &mut D)> {
        self.0.iter_mut()
    }

    /// returns the 'total degree' of a multidegree, i.e., the sum of the
    /// individual degree components
    pub fn total_deg(&self) -> D {
        self.0.values().sum()
    }

    /// returns the maximum index for which `self` contains an entry.
    ///
    /// in other words, this is the minimum value `n` for which we would
    /// consider `self` to be an element of the polynomial ring in `n` variables
    pub fn max_idx(&self) -> I {
        *self.0.keys().max().unwrap_or(&0)
    }
}

impl AddAssign<&Self> for MDeg {
    fn add_assign(&mut self, other: &Self) {
        for (&idx, &other_deg) in other.degs() {
            if let Some(self_deg_ref) = self.0.get_mut(&idx) {
                *self_deg_ref += other_deg;
            } else {
                self.0.insert(idx, other_deg);
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

impl fmt::Display for MDeg {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "MDeg{:?}", self.0)
    }
}

// term ------------------------------------------------------------------------

#[derive(Clone, PartialEq, Eq, Debug, Hash)]
pub struct Term<F: Field> {
    coef: F,
    mdeg: MDeg,
}

impl<F: Field> Term<F> {
    /// returns the term with the given coefficient and multidegree
    ///
    /// this is the most basic way one should create a new term struct. the
    /// current implementation is rather simply
    /// ```
    /// {
    ///    Term { coef, mdeg }
    /// }
    /// ```
    /// if, for some reason, term creation need be changed, this could help us
    /// ensure the change is uniform
    #[inline]
    pub fn from_coef_mdeg(coef: F, mdeg: MDeg) -> Self {
        Term { coef, mdeg }
    }

    /// returns a term with the given coefficient and multidegree `MDeg::0`
    /// 
    /// this is used for interpreting elements of the field as possible terms
    /// in polynomials over the field
    pub fn constant(coef: F) -> Self {
        Term::from_coef_mdeg(coef, MDeg::zero())
    }

    /// returns the constant term zero: `Term { coef: 0, mdeg: 0 }`  
    /// 
    /// this is a bit sketchy since the zero of the polynomial ring is
    /// typically not assigned a multidegree
    ///
    /// this would be sort of consistent with nonzero constant terms having
    /// multidegree zero, and if it doesn't cause any problems, then i'm
    /// willing to let the mathematics be a little uncomfortable
    pub fn zero() -> Self {
        Term::constant(F::ZERO)
    }

    /// returns the constant term one: `Term { coef: 1, mdeg: 0 }` 
    pub fn one() -> Self {
        Term::constant(F::ONE)
    }

    /// returns the monic term of given multidegree: `Term { coef: 1, mdeg }`
    pub fn monic(mdeg: MDeg) -> Self {
        Term::from_coef_mdeg(F::ONE, mdeg)
    }

    /// returns a term representing a single variable/indeterminate:
    /// - `var(0, k) = x_0^k = x^k`
    /// - `var(1, k) = x_1^k = y^k`
    /// - `var(2, k) = x_2^k = z^k`
    /// - `var(j, k) = x_j^k`
    pub fn var(idx: I, deg: D) -> Self {
        Term::monic(MDeg::single(idx, deg))
    }

    /// maps the polynomial element `self =: f ∈ F[x_1, ..., x_n]` to the
    /// corresponding
    /// polynomial function `eval_f: F^n -> F`, and the image of `x ∈ F^n`
    /// under this function is returned
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

/// the following multiplication operations are implemented:
///
/// term *= {&term, term}
/// term * term
/// &term * &term

impl<F: Field> MulAssign<&Self> for Term<F> {
    fn mul_assign(&mut self, other: &Self) {
        self.coef *= other.coef;
        self.mdeg += &other.mdeg;
    }
}

impl<F: Field> MulAssign for Term<F> {
    fn mul_assign(&mut self, other: Self) {
        *self *= &other;
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

impl<F: Field> Mul for &Term<F> {
    type Output = Term<F>;

    fn mul(self, other: Self) -> Self::Output {
        let mut out = self.clone();

        out *= other;

        // return
        out
    }
}

impl<F: Field> fmt::Display for Term<F> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}·X{:?}", self.coef, self.mdeg.0)
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
/// poly += {&term, term}
/// poly + {&term, term}
/// &poly + {&term, term}
/// poly += {&poly, poly}
/// poly + poly
/// &poly + &poly

impl<F: Field> AddAssign<&Term<F>> for Polynomial<F> {
    fn add_assign(&mut self, rhs: &Term<F>) {
        if let Some(term_ref) = self.terms.iter_mut().find(|t| t.mdeg == rhs.mdeg) {
            term_ref.coef += rhs.coef; // consider adding checks for zero coef
        } else {
            self.terms.push(rhs.clone());
        }
    }
}

impl<F: Field> AddAssign<Term<F>> for Polynomial<F> {
    fn add_assign(&mut self, rhs: Term<F>) {
        *self += &rhs;
    }
}

impl<F: Field> Add<&Term<F>> for Polynomial<F> {
    type Output = Self;

    fn add(mut self, rhs: &Term<F>) -> Self::Output {
        self += rhs;

        // return
        self
    }
}

impl<F: Field> Add<Term<F>> for Polynomial<F> {
    type Output = Self;

    fn add(self, rhs: Term<F>) -> Self::Output {
        self + &rhs
    }
}

impl<F: Field> Add<&Term<F>> for &Polynomial<F> {
    type Output = Polynomial<F>;

    fn add(self, rhs: &Term<F>) -> Self::Output {
        self.clone() + rhs
    }
}

impl<F: Field> Add<Term<F>> for &Polynomial<F> {
    type Output = Polynomial<F>;

    fn add(self, rhs: Term<F>) -> Self::Output {
        self + &rhs
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
/// poly *= {&term, term}
/// {&poly, poly} * {&term, term}
/// &poly * &poly
/// poly * poly

impl<F: Field> MulAssign<&Term<F>> for Polynomial<F> {
    fn mul_assign(&mut self, rhs: &Term<F>) {
        for term_ref in self.terms.iter_mut() {
            *term_ref *= rhs;
        }
    }
}

impl<F: Field> MulAssign<Term<F>> for Polynomial<F> {
    fn mul_assign(&mut self, rhs: Term<F>) {
        *self *= &rhs;
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

    fn mul(self, rhs: Term<F>) -> Self::Output {
        self * &rhs
    }
}

impl<F: Field> Mul<&Term<F>> for &Polynomial<F> {
    type Output = Polynomial<F>;

    fn mul(self, rhs: &Term<F>) -> Self::Output {
        self.clone() * rhs
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

impl<F: Field> fmt::Display for Polynomial<F> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut term_iter = self.terms.iter();

        if let Some(term) = term_iter.next() {
            write!(f, "{}", term)?;
        } else {
            return write!(f, "0_F[x]");
        }

        for term in term_iter {
            write!(f, " + {}", term)?;
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
        let zero = MDeg::zero();
        let ones = MDeg::ones(5);

        println!("zero = {}", zero);
        println!("ones = {}", ones);

        let a = MDeg::from_tuples(&[(0, 3), (3, 5), (4, 2)]);
        let b = MDeg::from_tuples(&[(0, 1), (1, 2), (4, 1), (5, 2)]);

        let c = MDeg::from_tuples(&[(0, 4), (1, 2), (3, 5), (4, 3), (5, 2)]);

        println!("\na = {}", a);
        println!("b = {}", b);

        let d = &a + &b;

        println!("\nc = {}", d);

        assert_eq!(a + b, c);

    }

    #[test]
    fn test_term() {
        let x = &Term::<f64>::var(0, 1);
        let y = &Term::<f64>::var(1, 1);
        let z = &Term::<f64>::var(2, 1);

        println!("x = {}", x);
        println!("y = {}", y);
        println!("z = {}", z);

        let t = &(x * y) * z;

        println!("\nt = x * y * z = {}", t);

        let c = Term::constant(37.0);

        println!("\nc = {}", c);

        let d = (&c * &c) * (&t * &t);

        println!("\nc^2 * t^2 = {}", d);
        println!("c^2 * t^2 = {:?}", d);

    }

    #[test]
    fn test_poly() {
        let c = |coef| Term::<f64>::constant(coef);
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

        let f = Polynomial::from_vec(vec![p, q, r]);

        println!("\nf = {}", f);

        let p = trm(5.0, [1, 2, 3]);
        let q = trm(7.0, [4, 0, 2]);
        let r = trm(13.0, [0, 5, 6]);

        let g = Polynomial::zero() + p + q + r;

        println!("\nf = {}", g);
    }
}