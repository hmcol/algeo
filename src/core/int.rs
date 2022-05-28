#![allow(unused)]

pub type Integer = i32;

pub fn gcd(mut a: Integer, mut b: Integer) -> Integer {
    let mut t: Integer;

    while b != 0 {
        t = b;
        b = a % b;
        a = t;
    }

    a.abs()
}

pub fn lcm(a: Integer, b: Integer) -> Integer {
    a / gcd(a, b) * b
}
