#![allow(unused)]


use crate::core::{
    element::{Element, Result},
    int::Integer,
};

use super::ring::{RingOps, Ring};

#[derive(Debug)]
pub enum Field {
    Q(Q),
    ZModPZ(ZModPZ),
}

pub trait FieldOps: RingOps {
    fn inv(&self, a: Element) -> Result<Element>;
    fn div(&self, a: Element, b: Element) -> Result<Element>;
}

#[derive(Debug)]
pub struct Q;

#[derive(Debug)]

pub struct ZModPZ(Integer);


pub trait FieldType: FieldOps {}

macro_rules! impl_ring_type {
    ($($t:ty)+) => {
        $(
            impl FieldType for $t {}
        )+
    };
}

impl_ring_type! { Q ZModPZ }

impl RingOps for Q {
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

impl FieldOps for Q {
    fn inv(&self, a: Element) -> Result<Element> {
        todo!()
    }

    fn div(&self, a: Element, b: Element) -> Result<Element> {
        todo!()
    }
}


impl RingOps for ZModPZ {
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

impl FieldOps for ZModPZ {
    fn inv(&self, a: Element) -> Result<Element> {
        todo!()
    }

    fn div(&self, a: Element, b: Element) -> Result<Element> {
        todo!()
    }
}