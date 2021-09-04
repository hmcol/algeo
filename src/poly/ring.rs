use std::fmt;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use xops::binop;

use crate::core::num::*;

// structs ---------------------------------------------------------------------

/// type for indexing the indeterminates
type I = usize;

/// type of the degrees in a multidegree
type D = i8;

/// Wrapper for constant values in the field `F`.
///
/// Maybe for implementing operations between all the types of elements in the polynomial ring `F[X]`, namely constants (here), terms, and polynomials.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct Const<F>(F);

/// wrapper for `Vec<P>` to represent multidegree of monomial terms in a
/// multivariate polynomial ring
#[derive(Clone, PartialEq, Eq, Default, Debug, Hash)]
pub struct MDeg(Vec<D>);

#[derive(Clone, PartialEq, Eq, Debug, Hash)]
pub struct Term<F: Field> {
    pub coef: F,
    pub mdeg: MDeg,
}

/// a polynomial with coefficients in the field F
#[derive(Clone, Debug, Hash)]
pub struct Polynomial<F: Field> {
    // pub? yeah pub.
    pub terms: Vec<Term<F>>,
}

/// possible way to handle interoperability between different types of elements in `F[X]`
enum _Ring<F: Field> {
    Const(F),
    Term(Term<F>),
    Poly(Polynomial<F>),
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
        }

        // return
        MDeg(vec)
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
        for deg in &self.0 {
            if *deg != 0 {
                return false;
            }
        }
        true
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
    pub fn from_pairs(pairs: &[(I, D)]) -> Self {
        /* let mut vec = Vec::new();
        let mut i = 0;

        for &(idx, deg) in pairs {
            if i == idx {
                vec.push(deg);
            } else {
                let diff = idx - i;
                vec.append(&mut vec![0; diff]);
                i += diff;
            }
        } */

        if let Some(n) = pairs.iter().map(|pair| pair.0).max() {
            let mut vec = vec![0; n + 1];

            for &(idx, deg) in pairs {
                vec[idx] = deg;
            }

            MDeg::from_vec(vec)
        } else {
            MDeg::zero()
        }
    }

    /// returns multidegree of all ones: `[1 ... 1]`
    ///
    /// used for the term `x_1` `x_2` \dots `x_n`, where each indeterminate
    /// occurs as a factor exactly once
    #[inline]
    pub fn ones(n: I) -> Self {
        MDeg::from_vec(vec![1; n])
    }

    /// returns the multidegree with the sole entry `idx: deg`
    ///
    /// can be used for representing variables/indeterminates:
    /// - the multidegree of `x = x_0` is taken to be `{0: 1}`
    /// - the multidegree of `y = x_1` is taken to be `{1: 1}`
    /// - the multidegree of `z = x_2` is taken to be `{2: 1}`
    #[inline]
    pub fn var(idx: I, deg: D) -> Self {
        MDeg::from_pairs(&[(idx, deg)])
    }

    /// returns iterator over the individual degree components
    ///
    /// Should replace return with in-house iterator struct.
    pub fn degs(&self) -> std::slice::Iter<D> {
        self.0.iter()
    }

    /// returns mutable iterator over the individual degree components
    ///
    /// Should replace return with in-house iterator struct.
    pub fn degs_mut(&mut self) -> std::slice::IterMut<D> {
        self.0.iter_mut()
    }

    /// returns the 'total degree' of a multidegree, i.e., the sum of the
    /// individual degree components
    pub fn total_deg(&self) -> D {
        self.0.iter().sum()
    }

    /// returns the maximum index for which `self` contains an entry.
    ///
    /// in other words, this is the minimum value `n` for which we would
    /// consider `self` to be an element of the polynomial ring in `n` variables
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
}

impl<F: Field> Term<F> {
    /// returns the term with the given coefficient and multidegree
    ///
    /// this is the most basic way one should create a new term struct. the
    /// current implementation is rather simply
    /// ```ignore
    /// {
    ///    Term { coef, mdeg }
    /// }
    /// ```
    /// if, for some reason, term creation need be changed, this could help us
    /// ensure the change is uniform
    #[inline]
    pub fn new(coef: F, mdeg: MDeg) -> Self {
        Term { coef, mdeg }
    }

    pub fn from_coef_degs(coef: F, degs: &[D]) -> Self {
        Term::new(coef, MDeg::from_vec(degs.into()))
    }

    /// returns a term with the given coefficient and multidegree `MDeg::0`
    ///
    /// this is used for interpreting elements of the field as possible terms
    /// in polynomials over the field
    pub fn constant(coef: F) -> Self {
        Term::new(coef, MDeg::zero())
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
        Term::new(F::ONE, mdeg)
    }

    /// returns a term representing a single variable/indeterminate:
    /// - `var(0, k) = x_0^k = x^k`
    /// - `var(1, k) = x_1^k = y^k`
    /// - `var(2, k) = x_2^k = z^k`
    /// - `var(j, k) = x_j^k`
    pub fn var(idx: I, deg: D) -> Self {
        Term::monic(MDeg::var(idx, deg))
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
        // should first check that self.coef is nonzero
        self.mdeg.is_succ(&other.mdeg)
    }
}

impl<F: Field> Polynomial<F> {
    pub fn from_vec(v: Vec<Term<F>>) -> Self {
        Polynomial { terms: v }
    }

    /// returns a polynomial with no terms
    pub fn zero() -> Self {
        Self::from_vec(Vec::new())
    }

    /// returns a polynomial with only the term `t`
    pub fn monomial(t: Term<F>) -> Self {
        Self::from_vec(vec![t])
    }

    /// returns the polynomial with only the constant term `1`
    pub fn one() -> Self {
        Self::monomial(Term::one())
    }

    /// number of (nonzero) terms
    pub fn term_count(&self) -> usize {
        self.terms.len()
    }

    /// Returns an iterator over (immutably) borrowed terms: `Item = &Term<F>`.
    ///
    /// Should replace return with in-house iterator struct.
    pub fn terms(&self) -> std::slice::Iter<Term<F>> {
        self.terms.iter()
    }

    /// Returns an iterator over mutably borrowed terms: `Item = &mut Term<F>`.
    ///
    /// Should replace return with in-house iterator struct.
    pub fn terms_mut(&mut self) -> std::slice::IterMut<Term<F>> {
        self.terms.iter_mut()
    }

    /// Evaluates `self` in `F[x_1, ..., x_n]` at the point `x` in `F^m`.
    ///
    /// # Panics
    ///
    /// Method panics if `m < n`, i.e., if `x` does not provide values up the
    /// maximum index variable of `self`.
    pub fn eval(&self, x: &[F]) -> F {
        self.terms().map(|t| t.eval(x)).sum()
    }
}

// defaults --------------------------------------------------------------------

impl<F: Field> Default for Polynomial<F> {
    fn default() -> Self {
        Self::zero()
    }
}

// conversions -----------------------------------------------------------------

impl<F: Field> From<F> for Term<F> {
    fn from(coef: F) -> Self {
        Term::constant(coef)
    }
}

impl<F: Field> From<Term<F>> for Polynomial<F> {
    fn from(term: Term<F>) -> Self {
        Polynomial::monomial(term)
    }
}

impl<F: Field> From<F> for Polynomial<F> {
    fn from(coef: F) -> Self {
        Polynomial::from(Term::from(coef))
    }
}

// operations ------------------------------------------------------------------

/// Given an implementation of `OpAssign<&B> for A`, this implements `OpAssign<B> for A` and `Op<B> for A`, with any specified arguments passed to `#[xops::binop(...)]`
///
/// this functionality will be moved to xops at some point
macro_rules! easy_ass {
    (@internal
        $TraitAss:ident $fn_ass:ident
        $Trait:ident $fn_:ident
        $(, $arg:ident)*
        for $Lhs:ty, $Rhs:ty
        $(where $($gens:tt)* )?
    ) => {
        impl$(< $($gens)* >)? $TraitAss<$Rhs> for $Lhs {
            fn $fn_ass(&mut self, rhs: $Rhs) {
                self.$fn_ass(&rhs);
            }
        }

        #[binop($($arg),*)]
        impl $(< $($gens)* >)? $Trait<$Rhs> for $Lhs {
            type Output = $Lhs;

            fn $fn_(mut self, rhs: $Rhs) -> Self::Output {
                self.$fn_ass(rhs);

                // return
                self
            }
        }
    };
    ($Lhs:ty [+ $(,$arg:ident)*] $($t:tt)*) => {
        easy_ass! { @internal AddAssign add_assign Add add $(,$arg)* for $Lhs, $($t)* }
    };
    ($Lhs:ty [- $(,$arg:ident)*] $($t:tt)*) => {
        easy_ass! { @internal SubAssign sub_assign Sub sub $(,$arg)* for $Lhs, $($t)* }
    };
    ($Lhs:ty [* $(,$arg:ident)*] $($t:tt)*) => {
        easy_ass! { @internal MulAssign mul_assign Mul mul $(,$arg)* for $Lhs, $($t)* }
    };
    ($Lhs:ty [/ $(,$arg:ident)*] $($t:tt)*) => {
        easy_ass! { @internal DivAssign div_assign Div div $(,$arg)* for $Lhs, $($t)* }
    };
}

impl AddAssign<&MDeg> for MDeg {
    fn add_assign(&mut self, rhs: &MDeg) {
        if self.len() < rhs.len() {
            let pad_zeros = &mut vec![0; std::cmp::max(0, rhs.len() - self.len())];
            self.0.append(pad_zeros);
        }

        self.degs_mut()
            .zip(rhs.degs())
            .for_each(|(deg_a, deg_b)| deg_a.add_assign(deg_b));
    }
}
easy_ass! { MDeg [+, refs_clone] MDeg }

impl SubAssign<&MDeg> for MDeg {
    fn sub_assign(&mut self, rhs: &MDeg) {
        if self.len() < rhs.len() {
            let pad_zeros = &mut vec![0; std::cmp::max(0, rhs.len() - self.len())];
            self.0.append(pad_zeros);
        }

        self.degs_mut()
            .zip(rhs.degs())
            .for_each(|(deg_a, deg_b)| deg_a.sub_assign(deg_b));
        
        self.trim_zeros();
    }
}
easy_ass! { MDeg [-, refs_clone] MDeg }

impl<F: Field> MulAssign<&Term<F>> for Term<F> {
    fn mul_assign(&mut self, rhs: &Term<F>) {
        self.coef *= rhs.coef;
        self.mdeg += &rhs.mdeg;
    }
}
easy_ass! { Term<F> [*, refs_clone] Term<F> where F: Field }

impl<F: Field> DivAssign<&Term<F>> for Term<F> {
    fn div_assign(&mut self, rhs: &Term<F>) {
        self.coef /= rhs.coef;
        self.mdeg -= &rhs.mdeg;
    }
}
easy_ass! { Term<F> [/, refs_clone] Term<F> where F: Field }

#[binop(derefs)]
impl<F: Field> Add for &Term<F> {
    type Output = Polynomial<F>;

    fn add(self, rhs: &Term<F>) -> Self::Output {
        Polynomial::from(self.clone()) + rhs
    }
}

impl<F: Field> AddAssign<&Term<F>> for Polynomial<F> {
    fn add_assign(&mut self, rhs: &Term<F>) {
        if let Some(i) = self.terms_mut().position(|t| t.mdeg == rhs.mdeg) {
            self.terms[i].coef += rhs.coef;

            // should be replaced with proper field zero checking
            // ya know, floating point error and stuff
            if self.terms[i].coef == F::ZERO {
                self.terms.remove(i);
            }
        } else {
            self.terms.push(rhs.clone());
        }
    }
}
easy_ass! { Polynomial<F> [+, refs_clone, commute] Term<F> where F: Field }

impl<F: Field> Neg for &Term<F> {
    type Output = Term<F>;

    fn neg(self) -> Self::Output {
        Term::new(-self.coef, self.mdeg.clone())
    }
}

impl<F: Field> SubAssign<&Term<F>> for Polynomial<F> {
    fn sub_assign(&mut self, rhs: &Term<F>) {
        if let Some(i) = self.terms_mut().position(|t| t.mdeg == rhs.mdeg) {
            self.terms[i].coef -= rhs.coef;

            // should be replaced with proper field zero checking
            // ya know, floating point error and stuff
            if self.terms[i].coef == F::ZERO {
                self.terms.remove(i);
            }
        } else {
            self.terms.push(-rhs);
        }
    }
}
easy_ass! { Polynomial<F> [-, refs_clone, commute] Term<F> where F: Field }

impl<F: Field> AddAssign<&Polynomial<F>> for Polynomial<F> {
    fn add_assign(&mut self, rhs: &Polynomial<F>) {
        for term in rhs.terms() {
            self.add_assign(term);
        }
    }
}
easy_ass! { Polynomial<F> [+, refs_clone] Polynomial<F> where F: Field }

impl<F: Field> SubAssign<&Polynomial<F>> for Polynomial<F> {
    fn sub_assign(&mut self, rhs: &Polynomial<F>) {
        for term in rhs.terms() {
            self.sub_assign(term);
        }
    }
}
easy_ass! { Polynomial<F> [-, refs_clone] Polynomial<F> where F: Field }

impl<F: Field> MulAssign<&Term<F>> for Polynomial<F> {
    fn mul_assign(&mut self, rhs: &Term<F>) {
        for term in self.terms_mut() {
            term.mul_assign(rhs);
        }
    }
}
easy_ass! { Polynomial<F> [*, refs_clone, commute] Term<F> where F: Field }

#[binop(derefs)]
impl<F: Field> Mul for &Polynomial<F> {
    type Output = Polynomial<F>;

    #[allow(clippy::suspicious_arithmetic_impl)]
    fn mul(self, other: &Polynomial<F>) -> Self::Output {
        other
            .terms
            .iter()
            .fold(Polynomial::zero(), |acc, t| acc + self * t)
    }
}

// shorthands ------------------------------------------------------------------

macro_rules! var_fn {
    (@with_doc $var:ident, $idx:literal, $doc_str:expr) => {
        #[doc = $doc_str]
        pub fn $var<F: Field>(deg: D) -> Term<F> {
            Term::var($idx, deg)
        }
    };
    (@doc_of $var:expr, $idx:expr) => {
        concat!(
            "Shorthand for `",
            $var,
            "`; i.e., the indeterminate of index `",
            $idx,
            "` in `F[X]`.\n\n",
            "As a Term, `",
            $var,
            "(d)` is monic with multidegree `{ ",
            $idx,
            ": d }`.",
            ""
        )
    };
    ($var:ident -> $idx:literal) => {
        var_fn! { @with_doc $var, $idx,
            var_fn!(@doc_of stringify!($var), stringify!($idx))
        }
    };
}

var_fn! { x -> 0 }
var_fn! { y -> 1 }
var_fn! { z -> 2 }
var_fn! { u -> 3 }
var_fn! { v -> 4 }
var_fn! { w -> 5 }

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
            }
        }
        Ok(())
    }
}

pub fn superscript(n: i8) -> String {
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
        n if n >= 10 => superscript(n / 10) + &superscript(n % 10),
        n if n < 0 => String::from("⁻") + &superscript(n.abs()),
        _ => String::from("[bad deg]")
    }
}

impl<F: Field> fmt::Display for Const<F> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.0)
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
            if term.coef < 0.0 {
                write!(f, " - {}", Term::new(-term.coef, term.mdeg.clone()))?;
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
        let zero = MDeg::zero();
        let ones = MDeg::ones(5);

        println!("zero = {}", zero);
        println!("ones = {}", ones);

        let a = MDeg::from_pairs(&[(0, 3), (3, 5), (4, 2)]);
        let b = MDeg::from_pairs(&[(0, 1), (1, 2), (4, 1), (5, 2)]);

        let c = MDeg::from_pairs(&[(0, 4), (1, 2), (3, 5), (4, 3), (5, 2)]);

        println!("\na = {}", a);
        println!("b = {}", b);

        let d = &a + &b;

        println!("\nc = {}", d);

        assert_eq!(a + b, c);
    }

    #[test]
    fn test_term() {
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
