use std::ops::{
    Add, AddAssign, Sub, SubAssign, Mul, MulAssign, Div, DivAssign, Index
};
use std::cmp::Ordering;
use std::fmt::{self, Write};

use super::num::*;



// multidegree -----------------------------------------------------------------

/// type of the powers/exponents/degrees in a multidegree
type P = i32;

/// wrapper for arrays to represent multidegree of monomial terms in a multivariate polynomial ring
/// 
/// e.g., monomial x0^2 x1^4 x2 x4^3 has multidegree (2, 4, 1, 0, 3)
#[derive(Copy, Clone, PartialEq, Eq, Debug, Hash)]
pub struct MDeg<const N: usize>(pub [P; N]);

impl<const N: usize> MDeg<N> {
    /// returns multidegree of all n's: `[n ... n]`
    const fn repeated(n: P) -> Self {
        MDeg([n; N])
    }

    /// returns multidegree of all zeros: `[0 ... 0]`
    /// 
    /// used for constant terms, i.e., terms without indeterminates
    pub const fn zeros() -> Self {
        Self::repeated(0)
    }

    /// returns multidegree of all ones: `[1 ... 1]`
    /// 
    /// used for the term `x_1` `x_2` \dots `x_N`, where each indeterminate occurs as a factor exactly once
    pub const fn ones() -> Self {
        Self::repeated(1)
    }

    /// returns iterator over the individual degree components
    pub fn degs(&self) -> std::slice::Iter<'_, P> {
        self.0.iter()
    }

    /// returns mutable iterator over the individual degree components
    pub fn degs_mut(&mut self) -> std::slice::IterMut<'_, P> {
        self.0.iter_mut()
    }

    /// returns the 'total degree' of a multidegree, i.e., the sum of the individual degree components
    pub fn total_deg(&self) -> P {
        self.degs().sum()
    }
}

/// implements the given operation & operation-assignment traits for
/// multidegrees componentwise
/// 
/// for example if the operation in question is star, `⋆`, then
/// `MDeg([a, b, c]) ⋆ MDeg([i, j, k])` gives `MDeg([a ⋆ i, b ⋆ j, c ⋆ k])
/// 
/// in particular, the operation-assignment trait is implemented first, with
/// the operation trait piggy-backing off it
macro_rules! impl_op_for_mdeg_componentwise {
    (
        $trait:ty,
        $func:ident,
        $ass_trait:ty,
        $ass_func:ident
    ) => {
        impl<const N: usize> $ass_trait for MDeg<N> {
            #[inline]
            fn $ass_func(&mut self, other: Self) {
                for (a_ref, b) in self.degs_mut().zip(other.degs()) {
                    a_ref.$ass_func(b);
                }
            }
        }

        impl<const N: usize> $trait for MDeg<N> {
            type Output = Self;
            
            #[inline]
            fn $func(mut self, other: Self) -> Self::Output {
                self.$ass_func(other);
        
                // return
                self
            }
        }
    };
}

impl_op_for_mdeg_componentwise! { Add, add, AddAssign, add_assign }
impl_op_for_mdeg_componentwise! { Sub, sub, SubAssign, sub_assign }


impl<const N: usize> Index<usize> for MDeg<N> {
    type Output = P;

    fn index(&self, index: usize) -> &Self::Output {
        self.0.index(index)
    }
}

impl<const N: usize> fmt::Display for MDeg<N> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut out = String::new();

        for d in self.degs() {
            write!(out, "{} ", d).unwrap();
        }
        out.pop();

        write!(f, "[{}]", out)
    }
}

// monomial orders -------------------------------------------------------------

/// enum for various possible monomial orders on multidegrees
pub enum MonomialOrder {
    Lex,
    RevLex,
    Grad,
    GrLex,
    GRevLex,
}

/// applies the given monomial order to a and b
pub fn cmp<const N: usize>(order: MonomialOrder, a: MDeg<N>, b: MDeg<N>) -> Ordering {
    match order {
        MonomialOrder::Lex => cmp_lex(a, b),
        MonomialOrder::RevLex => cmp_revlex(a, b),
        MonomialOrder::Grad => cmp_grad(a, b),
        MonomialOrder::GrLex => cmp_grlex(a, b),
        MonomialOrder::GRevLex => cmp_grevlex(a, b),
    }
}

/// lexicographic order
/// 
/// if a != b, compares first unequal degrees from the left
/// 
/// e.g., a < b iff ∃k s.t. a_k < b_k and a_i = b_i, for i = 0,...,k-1
pub fn cmp_lex<const N: usize>(a: MDeg<N>, b: MDeg<N>) -> Ordering {
    a.0.cmp(&b.0)
}

/// reverse lexicographic ordering on multidegrees
pub fn cmp_revlex<const N: usize>(a: MDeg<N>, b: MDeg<N>) -> Ordering {
    for (deg_a, deg_b) in a.degs().zip(b.degs()).rev() {
        match deg_a.cmp(deg_b) {
            Ordering::Equal => continue,
            lt_or_gt => return lt_or_gt,
        }
    }
    
    // return
    Ordering::Equal
}

/// graded order
/// 
/// simply compares the total degrees
/// 
/// i.e., the usual grading on a polynomial ring
pub fn cmp_grad<const N: usize>(a: MDeg<N>, b: MDeg<N>) -> Ordering {
    a.total_deg().cmp(&b.total_deg())
}

/// graded lexicographic order
/// 
/// applies the graded order; if equal, applies lexicographic
pub fn cmp_grlex<const N: usize>(a: MDeg<N>, b: MDeg<N>) -> Ordering {
    match cmp_grad(a, b) {
        Ordering::Equal => cmp_lex(a, b),
        lt_or_gt => lt_or_gt,
    }
}

/// graded reverse lexicographic order
/// 
/// applies the graded order; if equal, applies reverse lexicographic with the result negated
pub fn cmp_grevlex<const N: usize>(a: MDeg<N>, b: MDeg<N>) -> Ordering {
    match cmp_grad(a, b) {
        Ordering::Equal => cmp_revlex(a, b).reverse(),
        lt_or_gt => lt_or_gt,
    }
}


// term ------------------------------------------------------------------------

#[derive(Copy, Clone, PartialEq, Eq, Debug, Hash)]
pub struct Term<F: Field, const N: usize> {
    coef: F,
    mdeg: MDeg<N>,
}

impl<F: Field, const N: usize> Term<F, N> {
    pub fn from_coef_mdeg(coef: F, mdeg: MDeg<N>) -> Self {
        Term {
            coef,
            mdeg,
        }
    }

    /// returns a term with the given coefficient and multidegree `[0 ... 0]`
    pub fn constant(coef: F) -> Self {
        Term::from_coef_mdeg(coef, MDeg::zeros())
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
    pub fn monic(mdeg: MDeg<N>) -> Self {
        Term::from_coef_mdeg(F::ONE, mdeg)
    }

    /// maps the polynomial element `self` of `F[X]` to the corresponding polynomial function `F -> F`, given by evaluation, which is then evaluated at `x`
    pub fn eval(&self, x: &[F; N]) -> F {
        self.coef * x.iter()
            .zip(self.mdeg.degs())
            .map(|(x, d)| x.powi32(*d))
            .fold(F::ONE, |acc, x| acc * x) // should impl Product for Field
    }
}

/// implements the given operation & operation-assignment traits for
/// terms, with the coefficients being performed on by the very same operation
/// and the multidegrees being performed on by some corresponding operation
/// 
/// in particular, the operation-assignment trait is implemented first, with
/// the operation trait piggy-backing off it
macro_rules! impl_op_for_term {
    (
        $trait:ty,
        $func:ident,
        $ass_trait:ty,
        $ass_func:ident,
        $mdeg_ass_func:ident
    ) => {
        impl<F: Field, const N: usize> $ass_trait for Term<F, N> {
            fn $ass_func(&mut self, other: Self) {
                self.coef.$ass_func(other.coef);
                self.mdeg.$mdeg_ass_func(other.mdeg);
            }
        }

        impl<F: Field, const N: usize> $trait for Term<F, N> {
            type Output = Self;
        
            fn $func(mut self, other: Self) -> Self::Output {
                self.$ass_func(other);
        
                // return
                self
            }
        }
    };
}

impl_op_for_term! { Mul, mul, MulAssign, mul_assign, add_assign }
impl_op_for_term! { Div, div, DivAssign, div_assign, sub_assign }


impl<F: Field, const N: usize> fmt::Display for Term<F, N> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}X^{}", self.coef, self.mdeg)
    }
}

// polynomial ------------------------------------------------------------------

/// a polynomial over N variables with coefficients in the field F
#[derive(Clone, Debug, Hash)]
pub struct Polynomial<F:Field, const N: usize> {
    terms: Vec<Term<F, N>>,
}

impl<F: Field, const N: usize> Polynomial<F, N> {
    pub fn from_vec(v: Vec<Term<F, N>>) -> Self {
        Polynomial {
            terms: v,
        }
    }

    /// returns a polynomial with only the term `t`
    /// 
    /// not to be confused with the meaning of monomial which is the monic part of a term. not really sure what to do about terminology (get it?)
    pub fn monomial(t: Term<F, N>) -> Self {
        Self::from_vec(vec![t])
    }

    /// returns the 'soft-zero polynomial', i.e., the interpretation that the zero polynomial is 'the polynomial with no terms'
    /// 
    /// i might be more inclined to go with this one, since it is simple and should agree with a straightforward implementation of polynomial arithmetic.
    /// 
    /// for example if `p`, `q` in `F[X]`, then we might say their product pq is the sum over all terms u in p and v in q. if either p or q is the soft-zero polynomial, this will, again, produce the soft-zero polynomial
    /// 
    /// moreover it correctly does not assign `0` a multidegree
    pub fn zero() -> Self {
        Self::from_vec(Vec::new())
    }

    /// returns the 'hard-zero polynomial', i.e., the interpretation that the zero polynomial is 'the polynomial with only the zero term'
    /// 
    /// a possible argument in favor of this version would be that it actually incorporates the idea of a zero term, which may later help ground algorithms of polynomials with respect to terms.
    /// 
    /// idk what exactly, but maybe something with mutability where failure to act behaves more like `1` than `0`.
    /// 
    /// however, it does assign zero the multidegree `[0 ... 0]`
    pub fn zero_alt() -> Self {
        Self::monomial(Term::zero())
    }

    /// returns the polynomial with only the single term `1`, which correctly has the degree `[0 ... 0]`
    pub fn one() -> Self {
        Self::monomial(Term::one())
    }
    
    pub fn term_count(&self) -> usize {
        self.terms.len()
    }

    pub fn eval(&self, x: &[F; N]) -> F {
        self.terms
            .iter()
            .map(|t| t.eval(x))
            .fold(F::ZERO, |acc, y| acc + y) // should impl Sum for Field
    }
}


// uhhh idk which and how of the following to implement


impl<F: Field, const N: usize> AddAssign<Term<F,N>> for Polynomial<F, N> {
    fn add_assign(&mut self, rhs: Term<F,N>) {
        if let Some(term_ref) = self.terms
            .iter_mut()
            .find(|t| t.mdeg == rhs.mdeg)
        {
            term_ref.coef += rhs.coef;  // consider adding checks for zero coef
        } else {
            self.terms.push(rhs);
        }
    }
}

impl<F: Field, const N: usize> AddAssign for Polynomial<F, N> {
    fn add_assign(&mut self, other: Self) {
        for &t in &other.terms { *self += t };
    }
}

impl<F: Field, const N: usize> Add for Polynomial<F, N> {
    type Output = Self;

    fn add(mut self, other: Self) -> Self::Output {
        self += other;

        // return
        self
    }
}


impl<F: Field, const N: usize> Add for &Polynomial<F, N> {
    type Output = Polynomial<F, N>;

    fn add(self, other: Self) -> Self::Output {
        let mut out = self.clone();
        for &t in &other.terms { out += t };

        // return
        out
    }
}

impl<F: Field, const N: usize> Add<&Polynomial<F, N>> for Polynomial<F, N> {
    type Output = Polynomial<F, N>;

    fn add(self, other: &Self) -> Self::Output {
        &self + other
    }
}

impl<F: Field, const N: usize> Add<Polynomial<F, N>> for &Polynomial<F, N> {
    type Output = Polynomial<F, N>;

    fn add(self, other: Polynomial<F, N>) -> Self::Output {
        self + &other
    }
}







impl<F: Field, const N: usize> MulAssign<Term<F, N>> for Polynomial<F, N> {
    fn mul_assign(&mut self, other: Term<F, N>) {
        for term_ref in self.terms.iter_mut() {
            *term_ref *= other;
        }
    }
}

impl<F: Field, const N: usize> MulAssign<&Term<F, N>> for Polynomial<F, N> {
    fn mul_assign(&mut self, other: &Term<F, N>) {
        *self *= *other;
    }
}


impl<F: Field, const N: usize> Mul<Term<F, N>> for &Polynomial<F, N> {
    type Output = Polynomial<F, N>;

    fn mul(self, other: Term<F, N>) -> Self::Output {
        let mut out = self.clone();

        out *= other;

        // return
        out
    }
}

impl<F: Field, const N: usize> Mul<Term<F, N>> for Polynomial<F, N> {
    type Output = Self;

    fn mul(self, other: Term<F, N>) -> Self::Output {
        &self * other
    }
}

impl<F: Field, const N: usize> Mul<&Term<F, N>> for &Polynomial<F, N> {
    type Output = Polynomial<F, N>;

    fn mul(self, other: &Term<F, N>) -> Self::Output {
        self * *other
    }
}

impl<F: Field, const N: usize> Mul<&Term<F, N>> for Polynomial<F, N> {
    type Output = Polynomial<F, N>;

    fn mul(self, other: &Term<F, N>) -> Self::Output {
        self * *other
    }
}



impl<F: Field, const N: usize> Mul for &Polynomial<F, N> {
    type Output = Polynomial<F, N>;

    fn mul(self, other: Self) -> Self::Output {
        let mut out = Polynomial::zero();
        
        for t in &other.terms {
            out += self * t;
        }

        // return
        out
    }
}

impl<F: Field, const N: usize> Mul<&Polynomial<F, N>> for Polynomial<F, N> {
    type Output = Polynomial<F, N>;

    fn mul(self, other: &Self) -> Self::Output {
        &self * other
    }
}

impl<F: Field, const N: usize> Mul<Polynomial<F, N>> for &Polynomial<F, N> {
    type Output = Polynomial<F, N>;

    fn mul(self, other: Polynomial<F, N>) -> Self::Output {
        self * &other
    }
}

impl<F: Field, const N: usize> Mul for Polynomial<F, N> {
    type Output = Self;

    fn mul(self, other: Self) -> Self::Output {
        &self * &other
    }
}


impl<F: Field, const N: usize> MulAssign for Polynomial<F, N> {
    fn mul_assign(&mut self, other: Self) {
        *self = self.clone() * other
    }
}

impl<F: Field, const N: usize> MulAssign<&Self> for Polynomial<F, N> {
    fn mul_assign(&mut self, other: &Self) {
        *self = self.clone() * other
    }
}



impl<F: Field, const N: usize> fmt::Display for Polynomial<F, N> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut out = String::new();
        let mut term_iter = self.terms.iter();

        if let Some(t) = term_iter.next() {
            write!(out, "{}", t).unwrap();
        }

        for t in term_iter {
            //if t.coef != F::ZERO {
            write!(out, " + {}", t).unwrap();
            //} 
        }

        write!(f, "{}", out)
    }
}