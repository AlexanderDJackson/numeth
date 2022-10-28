/*****************************
 * Programming Assignment #2 *
 *****************************
 * Written by: Alex Jackson
 *
 * This project produces:
 *  - a lagrange polynomial for the points ((1, 2), (2, 1), (3, 3), (4, 2), (5, 3), (6, 4))
 *  - a lagrange polynomial for the function f(x) = e^(−x^2) * (2x + 1) using 6 chebyshev points
 *  - a cubic spine for the points ((1, 2), (2, 1), (3, 3), (4, 2), (5, 3), (6, 4))
 *
 * Note:
 *    Because we need to print the polynomials, I created a macro, which expands the points
 *    it takes as input into a lagrange polynomial. This means that before the program actually
 *    compiles, the compiler expands the macro, so at compile-time, there is no function call for
 *    the points, rather, an expression as a string.
 *
 *    For example, 
 *
 *      l!((1, 2), (2, 1), (3, 3))
 *
 *    would expand to:
 *
 *      2 * (x - 2) / (1 - 2) * (x - 3) / (1 - 3) + 1 * (x - 3) / (2 - 3) * (x - 1) / (2 - 1) + 3 * (x - 1) / (3 - 1) * (x - 2) / (3 - 2)
 *
 *    We then use another macro to evalutate this, as this is a string and I don't want to have an
 *    extra dependency to evaluate an expression if not necessary.
 *
 * Run Project:
 *      $ cargo run
 */

#[macro_use]
extern crate rulinalg;

use rulinalg::vector::Vector;
use rulinalg::matrix::{BaseMatrix, Matrix};

fn c(i: f64) -> f64 {
    2.0 * (((2.0 * i + 1.0) / 11.0) * std::f64::consts::PI).cos()
}

fn f(x: f64) -> f64 {
    std::f64::consts::E.powf(-x.powf(2.0)) * (x.powf(2.0) + 1.0)
}

fn dd(f: fn(f64) -> f64, a: &[f64]) -> f64 {
    match a.len() {
        0 => {
            std::f64::NAN
        },
        1 => {
            f(a[0])
        },
        2 => {
            (f(a[1]) - f(a[0])) / (a[1] - a[0])
        },
        _ => {
            (dd(f, &a[1..]) - dd(f, &a[0..a.len() - 1])) / (&a[a.len() - 1] - a[0])
        }
    }
}

fn pad(mut m: Matrix<f64>, i: u8, n: u8) -> Matrix<f64> {
    let n = n - m.cols() as u8 / 4;

    if i > 0 {
        let mut a: Matrix<f64> = matrix![0.0, 0.0, 0.0, 0.0];

        for _ in 0..i - 1 {
            a = a.hcat(&matrix![0.0, 0.0, 0.0, 0.0]);
        }

        a = a.hcat(&m);

        for _ in i..n - 1 {
            a = a.hcat(&matrix![0.0, 0.0, 0.0, 0.0]);
        }

        a
    } else {
        for _ in 0..n - 1 {
            m = m.hcat(&matrix![0.0, 0.0, 0.0, 0.0]); 
        }

        m
    }
}

fn spline(p: Vec<(f64, f64)>) -> Vector<f64> {
    let n = p.len() as u8;
    let mut a = pad(matrix![1.0, p[0].0, p[0].0.powf(2.0), p[0].0.powf(3.0)], 0, n);
    let mut b = vec![p[0].1];

    // Verify points are interpolated
    for i in 1..n - 1 {
        a = a.vcat(&pad(matrix![1.0, p[i as usize].0, p[i as usize].0.powf(2.0), p[i as usize].0.powf(3.0)], i, n));
        b.push(p[i as usize].1);
    }

    // Verify continuity
    for i in 0..n - 1 {
        a = a.vcat(&pad(matrix![1.0, p[i as usize + 1].0, p[i as usize + 1].0.powf(2.0), p[i as usize + 1].0.powf(3.0)], i, n));
        b.push(p[i as usize + 1].1);
    }

    // Verify differentiability
    for i in 0..n - 2 {
        a = a.vcat(&pad(matrix![0.0, 1.0, 2.0 * p[i as usize + 1].0, 3.0 * p[i as usize + 1].0.powf(2.0), 0.0, -1.0, -2.0 * p[i as usize + 1].0, -3.0 * p[i as usize + 1].0.powf(2.0)], i, n));
        b.push(0.0);
    }

    // Verify second differentiability
    for i in 0..n - 2 {
        a = a.vcat(&pad(matrix![0.0, 0.0, 2.0, 6.0 * p[i as usize + 1].0, 0.0, 0.0, -2.0, -6.0 * p[i as usize + 1].0], i, n));
        b.push(0.0);
    }

    // Verify second derivative of endpoints
    a = a.vcat(&pad(matrix![0.0, 0.0, 2.0, 6.0 * p[0].0], 0, n));
    b.push(0.0);

    a = a.vcat(&pad(matrix![0.0, 0.0, 2.0, 6.0 * p[n as usize - 1].0], n - 2, n));
    b.push(0.0);

    a.solve(Vector::new(b)).unwrap()
}

fn main() {
    let mut x;
    let mut fx;

    macro_rules! l {
        // One point
        // ((x, y))
        ((($x:expr, $y:expr))) => { $y };
        // x, y
        ($x:expr, $y:expr) => { l!((($x, $y))) };
        // (x, y)
        (($x:expr, $y:expr)) => { l!((($x, $y))) };

        // First term
        // ((x, y), (i, j), ...)
        // We return the body and separate the (x, y) to denote that we've already used it
        ((($x:expr, $y:expr), $(($i:expr, $j:expr)),+)) => {
             $y $(* (x - $i) / ($x - $i))+ + l!(($(($i, $j)),+); (($x, $y)))
        };

        // Nth term
        // ((x, y), (i, j)); ((a, b))
        // We always use the first term as the concerning point, then move it to the end so we don't reuse it
        ((($x:expr, $y:expr), $(($i:expr, $j:expr)),+); ($(($a:expr, $b:expr)),+)) => {
         $y $(* (x - $i) / ($x - $i))+ $(* (x - $a) / ($x - $a))+ + l!(($(($i, $j)),+); ($(($a, $b)),+, ($x, $y)))
        };

        // Last term
        // ((x, y)); ((a, b), ...)
        ((($x:expr, $y:expr)); ($(($a:expr, $b:expr)),+)) => {
            $y $(* (x - $a) / ($x - $a))+ 
        };
    }

    /*
    macro_rules! n {
        // One point
        // ((x, y))
        ((($x:expr, $y:expr))) => { $y };
        // x, y
        ($x:expr, $y:expr) => { n!((($x, $y))) };
        // (x, y)
        (($x:expr, $y:expr)) => { n!((($x, $y))) };

        // First term
        // ((x, y), (i, j), ...)
        // We return the body and separate the (x, y) to denote that we've already used it
        ((($x:expr, $y:expr), $(($i:expr, $j:expr)),+)) => {
             $y + n!(((($x, $y); $(($i, $j)),+))
        };

        // Nth term
        // ((x, y), (i, j)); ((a, b))
        // We always use the first term as the concerning point, then move it to the end so we don't reuse it
        (($(($x:expr, $y:expr)),+); ($(($a:expr, $b:expr)),+)) => {
         $y $(* (x - $i) / ($x - $i))+ $(* (x - $a) / ($x - $a))+ + l!(($(($i, $j)),+); ($(($a, $b)),+, ($x, $y)))
        };

        // Last term
        // ((x, y)); ((a, b), ...)
        ((($x:expr, $y:expr)); ($(($a:expr, $b:expr)),+)) => {
            $y $(* (x - $a) / ($x - $a))+ 
        };
    }
    */

    println!("Lagrange polynomial for points: (1.0, 2.0), (2.0, 1.0), (3.0, 3.0), (4.0, 2.0), (5.0, 3.0), (6.0, 4.0)");

    // Evaluating the lagrange polynomial for part 1
    for i in [1.0, 1.5, 2.25, 3.25, 4.621].iter() {
        x = i;

        fx = l!(
            (
                (1.0, 2.0),
                (2.0, 1.0),
                (3.0, 3.0),
                (4.0, 2.0),
                (5.0, 3.0),
                (6.0, 4.0)
            )
        );

        println!("P({}) = {}", x, fx);
    }

    println!("\nChebyshev polynomial for e^-x^2 * (x^2 + 1)");

    // Evaluating the chebyshev polynomial for part 2
    for i in [1.0, 1.5, 2.25, 3.25, 4.621].iter() {
        x = i;

        fx = l!(
            (
                (c(0.0), f(c(0.0))),
                (c(1.0), f(c(1.0))),
                (c(2.0), f(c(2.0))),
                (c(3.0), f(c(3.0))),
                (c(4.0), f(c(4.0))),
                (c(5.0), f(c(5.0))),
                (c(6.0), f(c(6.0)))
            )
        );

        println!("P({}) = {}", x, fx);
    }

    println!("{}", spline(vec![(1.0, 2.0), (2.0, 1.0), (3.0, 3.0), (4.0, 2.0), (5.0, 3.0), (6.0, 4.0)]));
}
