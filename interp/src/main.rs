/*****************************
 * Programming Assignment #2 *
 *****************************
 * Written by: Alex Jackson
 *
 * This project produces:
 *  - a lagrange polynomial for the points ((1, 2), (2, 1), (3, 3), (4, 2), (5, 3), (6, 4))
 *  - a lagrange polynomial for the function f(x) = e^(−x^2) * (2x + 1) using 6 chebyshev nodes
 *  - a cubic spine for the points ((1, 2), (2, 1), (3, 3), (4, 2), (5, 3), (6, 4))
 *
 * Note:
 *    Because I thought we needed to print out the polynomials, I created a macro, which expands the points
 *    it takes as input into a lagrange polynomial. This means that before the program actually
 *    compiles, the compiler expands the macro, so at compile-time, there is no function call for
 *    the points, but rather an expression.
 *
 *    For example, 
 *
 *      l!((1, 2), (2, 1), (3, 3))
 *
 *    would expand to:
 *
 *      2 * (x - 2) / (1 - 2) * (x - 3) / (1 - 3) + 1 * (x - 3) / (2 - 3) * (x - 1) / (2 - 1) + 3 * (x - 1) / (3 - 1) * (x - 2) / (3 - 2)
 *
 * Run Project:
 *      $ cargo run
 */

#[macro_use]
extern crate rulinalg;

use rulinalg::vector::Vector;
use rulinalg::matrix::{BaseMatrix, Matrix};

// Given i, this function returns the ith chebyshev node
// Note: i must be in [1, n], otherwise polynomial interpolation could fail
fn c(i: f64) -> f64 {
    2.0 * (((2.0 * i - 1.0) / 12.0) * std::f64::consts::PI).cos()
}

// Given x, this function returns the e^-x^2 * (x^2 + 1)
fn f(x: f64) -> f64 {
    std::f64::consts::E.powf(-x.powf(2.0)) * (x.powf(2.0) + 1.0)
}

// Given a function and a set of x values, this function returns the
// divided diffence of the x values: f[x0..xn]
/*
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
*/

/* Given:
 *  - m: contains a constraint for some spline function
 *  - i: specifies the spline function
 *  - n: specifies the total number of spline functions
 * Return:
 *  - A matrix that pads the constraint with zeroes to 
 *    associate it with the correct spline function
 */
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

// Given a vector of points, sorted by x value, to interpolate, this function returns
// a vector of coefficients for functions to build a piecewise natural cubic spline function
fn spline_solve(p: &Vec<(f64, f64)>) -> Vector<f64> {
    let n = p.len() as u8;
    let mut a = pad(matrix![1.0, p[0].0, p[0].0.powf(2.0), p[0].0.powf(3.0)], 0, n);
    let mut b = vec![p[0].1];

    // Points interpolation constraint
    for i in 1..n - 1 {
        a = a.vcat(&pad(matrix![1.0, p[i as usize].0, p[i as usize].0.powf(2.0), p[i as usize].0.powf(3.0)], i, n));
        b.push(p[i as usize].1);
    }

    // Continuity constraint
    for i in 0..n - 1 {
        a = a.vcat(&pad(matrix![1.0, p[i as usize + 1].0, p[i as usize + 1].0.powf(2.0), p[i as usize + 1].0.powf(3.0)], i, n));
        b.push(p[i as usize + 1].1);
    }

    // Differentiability constraint
    for i in 0..n - 2 {
        a = a.vcat(&pad(matrix![0.0, 1.0, 2.0 * p[i as usize + 1].0, 3.0 * p[i as usize + 1].0.powf(2.0), 0.0, -1.0, -2.0 * p[i as usize + 1].0, -3.0 * p[i as usize + 1].0.powf(2.0)], i, n));
        b.push(0.0);
    }

    // Second differentiability constraint
    for i in 0..n - 2 {
        a = a.vcat(&pad(matrix![0.0, 0.0, 2.0, 6.0 * p[i as usize + 1].0, 0.0, 0.0, -2.0, -6.0 * p[i as usize + 1].0], i, n));
        b.push(0.0);
    }

    // Second derivative endpoints constraint
    a = a.vcat(&pad(matrix![0.0, 0.0, 2.0, 6.0 * p[0].0], 0, n));
    b.push(0.0);

    a = a.vcat(&pad(matrix![0.0, 0.0, 2.0, 6.0 * p[n as usize - 1].0], n - 2, n));
    b.push(0.0);

    // Return the vector containing values of a, b, c, and d
    // that correspond to the appropriate spline function
    a.solve(Vector::new(b)).unwrap()
}

/* Given:
 *  - x: the value at which to evaluate the function
 *  - v: the coefficients of the spline functions
 *  - p: the points the splines interpolate
 * Returns:
 *  - the value of the natural cubic spline at x
 */
fn spline_evaluate(x: f64, v: &Vec<f64>, p: &Vec<(f64, f64)>) -> f64 {
    // The range in which x lies
    // these are the indexes of the point, so p[0], not actually 0
    let mut a = 0;
    let mut b = 0;

    for (i, point) in p.iter().enumerate() {
        if point.0 > p[a].0 && point.0 < x {
            a = i;
        } else if point.0 < p[b].0 && point.0 > x {
            b = 1
        }
    }

    let mut y: f64 = 0.0;

    for i in (a * 4)..(a * 4 + 4) {
        y += v[i] * x.powf((i - (a * 4)) as f64);
    }

    y
}

fn main() {
    let mut x;
    let mut fx;
    let points = [1.0, 1.5, 2.25, 3.25, 4.621];

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

    println!("Lagrange polynomial for points: (1.0, 2.0), (2.0, 1.0), (3.0, 3.0), (4.0, 2.0), (5.0, 3.0), (6.0, 4.0)");

    // Evaluating the lagrange polynomial for part 1
    for i in points {
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

    println!("\nLagrange polynomial for e^-x^2 * (x^2 + 1) using six Chebyshev nodes");

    // Evaluating the chebyshev polynomial for part 2
    for i in points {
        x = i;

        fx = l!(
            (
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

    println!("\nNatural Cubic Spline for (1.0, 2.0), (2.0, 1.0), (3.0, 3.0), (4.0, 2.0), (5.0, 3.0), (6.0, 4.0)");

    let p = vec![(1.0, 2.0), (2.0, 1.0), (3.0, 3.0), (4.0, 2.0), (5.0, 3.0), (6.0, 4.0)];
    let coefficients = spline_solve(&p).into_vec();

    // print the cubic spline
    for (i, c) in coefficients.iter().enumerate() {
        match i % 4 {
            0 => { if i / 4 == (p.len() - 1) / 2 { print!("S(x) =") }; print!("\tS{}(x) = {}", i / 4, c); },
            1 => { if c >= &0.0 { print!(" + {}x", c); } else { print!(" - {}x", -1.0 * c); } },
            2 => { if c >= &0.0 { print!(" + {}x^2", c); } else { print!(" - {}x^2", -1.0 * c); } },
            3 => { if c >= &0.0 { println!(" + {}x^3", c); } else { println!(" - {}x^3", -1.0 * c); } },
            _ => { todo!() }
        }
    }

    println!();

    // print the spline evaluated at the required points
    for i in points {
        println!("S({}) = {}", i, spline_evaluate(i, &coefficients, &p));
    }
}
