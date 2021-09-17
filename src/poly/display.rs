use std::fmt::{Display, Formatter, Result};

use super::{
    mdeg::MultiDegree,
    elts::{Polynomial, Term},
};
use crate::core::num::Field;

// helpers ---------------------------------------------------------------------

impl MultiDegree {
    fn write_var(&self, f: &mut Formatter, idx: usize, ident: &str) -> Result {
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

// implementations -------------------------------------------------------------

impl Display for MultiDegree {
    fn fmt(&self, f: &mut Formatter) -> Result {
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

impl<F: Field> Display for Term<F> {
    fn fmt(&self, f: &mut Formatter) -> Result {
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

impl Display for Polynomial<f64> {
    fn fmt(&self, f: &mut Formatter) -> Result {
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
        Result::Ok(())
    }
}

impl Display for Polynomial<crate::core::frac::Frac> {
    fn fmt(&self, f: &mut Formatter) -> Result {
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
        Result::Ok(())
    }
}
