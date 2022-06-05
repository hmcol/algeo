use super::{int::Integer, matrix::Mat, rat::Rational};

/// main type for representing mathematical elements of domains

#[derive(Debug, Clone)]
pub enum Element {
    Int(Int),
    IntPair(IntPair),
    IntVec(IntVec),
    Rat(Rat),
    Letter(Letter),
    Word(Word),
    ElementPair(ElementPair),
    ElementVec(ElementVec),
    ElementMat(ElementMat),
}

#[allow(non_snake_case)]
pub fn Elt<D: ElementDataType>(data: D) -> Element {
    data.into_element()
}

#[allow(non_snake_case)]
pub fn Pair(a: Element, b: Element) -> ElementPair {
    ElementPair::new(a, b)
}

/// most basic mathematical integer type
///
/// any rust integer which is meant to be interpreted mathematically must go
/// through this type
///
/// the integer operations here are not to be confused with those defined for
/// the various structures which go by the same name (e.g., the abelian group Z
/// and the ring Z)
#[derive(Debug, Clone)]
pub struct Int(pub Integer);

#[derive(Debug, Clone)]
pub struct IntPair(pub Integer, pub Integer);

#[derive(Debug, Clone)]
pub struct IntVec(pub Vec<Integer>);

#[derive(Debug, Clone)]
pub struct Rat(pub Rational);

#[derive(Debug, Clone)]
pub struct Letter(char);

#[derive(Debug, Clone)]
pub struct Word(String);

#[derive(Debug, Clone)]
pub struct ElementPair(Box<Element>, Box<Element>);

impl ElementPair {
    pub fn new(left: Element, right: Element) -> Self {
        ElementPair(Box::new(left), Box::new(right))
    }
    pub fn into_tuple(self) -> (Element, Element) {
        (*self.0, *self.1)
    }
}

#[derive(Debug, Clone)]
pub struct ElementVec(Vec<Element>);

#[derive(Debug, Clone)]
pub struct ElementMat(pub Mat<Element>);

#[derive(Debug, Default)]
pub struct DataTypeError(String);

pub type Result<T> = std::result::Result<T, DataTypeError>;
pub trait ElementDataType: Sized {
    fn read_data(e: Element) -> Result<Self>;
    fn into_element(self) -> Element;
}

macro_rules! impl_edt {
    ($t:ident) => {
        impl ElementDataType for $t {
            fn read_data(e: Element) -> Result<Self> {
                match e {
                    Element::$t(data) => Ok(data),
                    e => Err(DataTypeError(format!(
                        "no! something expected `Element::{}` but found `Element::{e:?}`",
                        stringify!($t)
                    ))),
                }
            }
            fn into_element(self) -> Element {
                Element::$t(self)
            }
        }
    };
    ($t:ident $($s:ident)+) => {
        impl_edt! { $t }
        impl_edt! { $($s)+ }
    };
}

impl_edt! { Int IntPair IntVec Rat Letter Word ElementPair ElementVec ElementMat }

#[cfg(test)]
mod tests {
    use crate::core::element::{ElementDataType, Word};

    use super::{Element, Int};

    #[test]
    fn elt() {
        println!();

        let elt_a = Element::Int(Int(34544));
        println!("elt_a = {elt_a:?}");
        let read_a = Int::read_data(elt_a);
        println!("read_a = {read_a:?}\n");

        let elt_b = Element::Word(Word("syzygy".to_string()));
        println!("elt_b = {elt_b:?}");
        let read_b = Int::read_data(elt_b);
        println!("read_b = {read_b:?}\n");

        let elt_b2 = Element::Word(Word("syzygy".to_string()));
        let read_b2 = Word::read_data(elt_b2);
        println!("read_b2 = {read_b2:?}\n");
    }
}
