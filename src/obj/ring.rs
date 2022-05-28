#![allow(unused)]

use crate::core::{element::*, matrix::Mat};

/// object in the category of rings
#[derive(Debug)]
pub enum Ring {
    Z(Z),
    ZModNZ(ZModNZ),
    Pair(PairRing),
    Prod(ProdRing),
    Matrix(MatrixRing),
}

pub trait RingOps {
    fn zero(&self) -> Element;
    fn add(&self, a: Element, b: Element) -> Result<Element>;
    fn sub(&self, a: Element, b: Element) -> Result<Element>;
    fn neg(&self, a: Element) -> Result<Element>;
    fn one(&self) -> Element;
    fn mul(&self, a: Element, b: Element) -> Result<Element>;
}

impl RingOps for Ring {
    fn zero(&self) -> Element {
        match self {
            Ring::Z(r) => r.zero(),
            Ring::ZModNZ(r) => r.zero(),
            Ring::Pair(r) => r.zero(),
            Ring::Prod(_) => todo!(),
            Ring::Matrix(_) => todo!(),
        }
    }

    fn add(&self, a: Element, b: Element) -> Result<Element> {
        match self {
            Ring::Z(z) => z.add(a, b),
            Ring::ZModNZ(_) => todo!(),
            Ring::Pair(rxs) => rxs.add(a, b),
            Ring::Prod(_) => todo!(),
            Ring::Matrix(mat_r) => mat_r.add(a, b),
        }
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

trait RingType: RingOps {}

macro_rules! impl_ring_type {
    ($($t:ty)+) => {
        $(
            impl RingType for $t {}
        )+
    };
}

impl_ring_type! { Z ZModNZ PairRing ProdRing MatrixRing }

#[derive(Debug)]
pub struct Z;

#[derive(Debug)]
pub struct ZModNZ(i32);

#[derive(Debug)]
pub struct PairRing(Box<Ring>, Box<Ring>);

#[derive(Debug)]
pub struct ProdRing(Vec<Ring>);

#[derive(Debug)]
pub struct MatrixRing {
    coef: Box<Ring>,
    rows: usize,
    cols: usize,
}



impl RingOps for Z {
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

impl RingOps for ZModNZ {
    fn zero(&self) -> Element {
        todo!()
    }

    fn add(&self, a: Element, b: Element) -> Result<Element> {
        todo!()
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

impl RingOps for PairRing {
    fn zero(&self) -> Element {
        Elt(ElementPair::new(self.0.zero(), self.1.zero()))
    }

    fn add(&self, a: Element, b: Element) -> Result<Element> {
        let a = ElementPair::read_data(a)?.into_tuple();
        let b = ElementPair::read_data(b)?.into_tuple();

        Ok(Elt(ElementPair::new(
            self.0.add(a.0, b.0)?,
            self.1.add(a.1, b.1)?,
        )))
    }

    fn sub(&self, a: Element, b: Element) -> Result<Element> {
        todo!()
    }

    fn neg(&self, a: Element) -> Result<Element> {
        todo!()
    }

    fn one(&self) -> Element {
        Elt(ElementPair::new(self.0.one(), self.1.one()))
    }

    fn mul(&self, a: Element, b: Element) -> Result<Element> {
        todo!()
    }
}

impl RingOps for ProdRing {
    fn zero(&self) -> Element {
        todo!()
    }

    fn add(&self, a: Element, b: Element) -> Result<Element> {
        todo!()
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

impl RingOps for MatrixRing {
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

