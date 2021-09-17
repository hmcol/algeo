#![allow(unused)]

/// Implements the `Default` trait by calling `Self::zero`.
///
/// For hopefully obvious reasons, this only works if the given type has an unambiguous function called `zero`.
// #[macro_export]
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

macro_rules! fn_var_term {
    (@with_doc $var:ident, $idx:literal, $doc_str:expr) => {
        #[doc = $doc_str]
        pub fn $var<F: Field>(deg: u8) -> Term<F> {
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

macro_rules! fn_vars {
    (@idx x) => { 0 };
    (@idx y) => { 1 };
    (@idx z) => { 2 };
    (@idx w) => { 3 };
    (@idx u) => { 4 };
    (@idx v) => { 5 };
    ($t:ty: $($var:ident)*) => {
        $(
            fn $var(deg: u8) -> Polynomial<$t> {
                crate::poly::elts::Polynomial::<$t>::var(fn_vars!(@idx $var), deg)
            }
        )*
    };
}
