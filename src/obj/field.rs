#![allow(unused)]

use crate::core::{
    element::{DataTypeError, Element, ElementDataType, Elt, Rat, Result},
    int::Integer,
};

use super::{ring::Ring, Q};

pub trait Field: Ring {
    fn inv(&self, a: Element) -> Result<Element>;
    fn div(&self, a: Element, b: Element) -> Result<Element>;
}

pub type FieldBox = Box<dyn Field>;

impl Field for Q {
    fn inv(&self, a: Element) -> Result<Element> {
        let a = Rat::read_data(a)?.0;

        match a.recip() {
            Some(a_inv) => Ok(Elt(Rat(a_inv))),
            // idk if this is the right sort of error
            // should be like "divide by zero" error
            None => Err(DataTypeError::default()),
        }
    }

    fn div(&self, a: Element, b: Element) -> Result<Element> {
        let a = Rat::read_data(a)?.0;
        let b = Rat::read_data(b)?.0;

        match b.recip() {
            Some(b_inv) => Ok(Elt(Rat(a * b_inv))),
            // idk if this is the right sort of error
            // should be like "divide by zero" error
            None => Err(DataTypeError::default()),
        }
    }
}

pub struct FieldProd(Vec<FieldBox>);
