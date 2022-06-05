#![allow(unused)]

use super::{ZModNZ, Q, Z};
use crate::core::{element::*, matrix::Mat, rat::Rational};

pub trait Ring {
    fn zero(&self) -> Element;
    fn add(&self, a: Element, b: Element) -> Result<Element>;
    fn sub(&self, a: Element, b: Element) -> Result<Element>;
    fn neg(&self, a: Element) -> Result<Element>;
    fn one(&self) -> Element;
    fn mul(&self, a: Element, b: Element) -> Result<Element>;
}

pub type RingBox = Box<dyn Ring>;

impl Ring for Z {
    fn zero(&self) -> Element {
        Elt(Int(0))
    }

    fn add(&self, a: Element, b: Element) -> Result<Element> {
        let a = Int::read_data(a)?.0;
        let b = Int::read_data(b)?.0;

        Ok(Elt(Int(a + b)))
    }

    fn sub(&self, a: Element, b: Element) -> Result<Element> {
        let a = Int::read_data(a)?.0;
        let b = Int::read_data(b)?.0;

        Ok(Elt(Int(a - b)))
    }

    fn neg(&self, a: Element) -> Result<Element> {
        let a = Int::read_data(a)?.0;

        Ok(Elt(Int(-a)))
    }

    fn one(&self) -> Element {
        Elt(Int(1))
    }

    fn mul(&self, a: Element, b: Element) -> Result<Element> {
        let a = Int::read_data(a)?.0;
        let b = Int::read_data(b)?.0;

        Ok(Elt(Int(a * b)))
    }
}

impl Ring for Q {
    fn zero(&self) -> Element {
        Elt(Rat(Rational::ZERO))
    }

    fn add(&self, a: Element, b: Element) -> Result<Element> {
        let a = Rat::read_data(a)?.0;
        let b = Rat::read_data(b)?.0;

        Ok(Elt(Rat(a + b)))
    }

    fn sub(&self, a: Element, b: Element) -> Result<Element> {
        let a = Rat::read_data(a)?.0;
        let b = Rat::read_data(b)?.0;

        Ok(Elt(Rat(a - b)))
    }

    fn neg(&self, a: Element) -> Result<Element> {
        let a = Rat::read_data(a)?.0;

        Ok(Elt(Rat(-a)))
    }

    fn one(&self) -> Element {
        Elt(Rat(Rational::ONE))
    }

    fn mul(&self, a: Element, b: Element) -> Result<Element> {
        let a = Rat::read_data(a)?.0;
        let b = Rat::read_data(b)?.0;

        Ok(Elt(Rat(a * b)))
    }
}

impl Ring for ZModNZ {
    fn zero(&self) -> Element {
        Elt(Int(0))
    }

    fn add(&self, a: Element, b: Element) -> Result<Element> {
        let a = Int::read_data(a)?.0;
        let b = Int::read_data(b)?.0;

        let out = (a + b) % self.modulus;

        Ok(Elt(Int(out)))
    }

    fn sub(&self, a: Element, b: Element) -> Result<Element> {
        let a = Int::read_data(a)?.0;
        let b = Int::read_data(b)?.0;

        let out = (a - b) % self.modulus;

        Ok(Elt(Int(out)))
    }

    fn neg(&self, a: Element) -> Result<Element> {
        let a = Int::read_data(a)?.0;

        let out = (-a) % self.modulus;

        Ok(Elt(Int(out)))
    }

    fn one(&self) -> Element {
        Elt(Int(1))
    }

    fn mul(&self, a: Element, b: Element) -> Result<Element> {
        let a = Int::read_data(a)?.0;
        let b = Int::read_data(b)?.0;

        let out = (a * b) % self.modulus;

        Ok(Elt(Int(out)))
    }
}

pub struct RingProd(RingBox, RingBox);

impl Ring for RingProd {
    fn zero(&self) -> Element {
        Elt(Pair(self.0.zero(), self.1.zero()))
    }

    fn add(&self, a: Element, b: Element) -> Result<Element> {
        let a = ElementPair::read_data(a)?.into_tuple();
        let b = ElementPair::read_data(b)?.into_tuple();

        let add_0 = self.0.add(a.0, b.0)?;
        let add_1 = self.1.add(a.1, b.1)?;

        Ok(Elt(Pair(add_0, add_1)))
    }

    fn sub(&self, a: Element, b: Element) -> Result<Element> {
        let a = ElementPair::read_data(a)?.into_tuple();
        let b = ElementPair::read_data(b)?.into_tuple();

        let sub_0 = self.0.sub(a.0, b.0)?;
        let sub_1 = self.1.sub(a.1, b.1)?;

        Ok(Elt(Pair(sub_0, sub_1)))
    }

    fn neg(&self, a: Element) -> Result<Element> {
        let a = ElementPair::read_data(a)?.into_tuple();

        let neg_0 = self.0.neg(a.0)?;
        let neg_1 = self.1.neg(a.1)?;

        Ok(Elt(Pair(neg_0, neg_1)))
    }

    fn one(&self) -> Element {
        Elt(Pair(self.0.one(), self.1.one()))
    }

    fn mul(&self, a: Element, b: Element) -> Result<Element> {
        let a = ElementPair::read_data(a)?.into_tuple();
        let b = ElementPair::read_data(b)?.into_tuple();

        let mul_0 = self.0.mul(a.0, b.0)?;
        let mul_1 = self.1.mul(a.1, b.1)?;

        Ok(Elt(Pair(mul_0, mul_1)))
    }
}

/// is this worth putting here? idk
pub struct MatrixRing {
    coef: Box<dyn Ring>,
    rows: usize,
    cols: usize,
}

impl Ring for MatrixRing {
    fn zero(&self) -> Element {
        Elt(ElementMat(Mat::new(
            self.rows,
            self.cols,
            vec![self.coef.zero(); self.rows * self.cols],
        )))
    }

    fn add(&self, a: Element, b: Element) -> Result<Element> {
        let mat_a = ElementMat::read_data(a)?.0;
        let mat_b = ElementMat::read_data(b)?.0;

        if mat_a.size() != mat_b.size() {
            // idk if this should be treated the same
            return Err(DataTypeError::default());
        }

        let mut out = Vec::new();

        for (a_ij, b_ij) in mat_a.entries().iter().zip(mat_b.entries().iter()) {
            // reason we can't just do iterators here is that this `?` doesn't
            // work inside a closure
            out.push(self.coef.add(a_ij.clone(), b_ij.clone())?);
        }

        let m = Mat::new(mat_a.rows(), mat_a.cols(), out);

        Ok(Elt(ElementMat(m)))
    }

    fn sub(&self, a: Element, b: Element) -> Result<Element> {
        todo!()
    }

    fn neg(&self, a: Element) -> Result<Element> {
        todo!()
    }

    fn one(&self) -> Element {
        todo!()
    }

    fn mul(&self, a: Element, b: Element) -> Result<Element> {
        todo!()
    }
}
