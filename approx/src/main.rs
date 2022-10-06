/*
 * Written By: Alex Jackson (asj18a@acu.edu)
 *
 * Implementation of Approximation Methods:
 *  - Bisection method
 *  - Newton’s method
 *  - Secant method
 *  - Muller’s method
 *  - Inverse quadratic iteration
 *  - False position method (Bracketed secant method)
 *
 * Running the program:
 *  - Ensure you have the latest rust toolchain
 *  - Run the following command: cargo run
 *
 *  */

#[macro_use]
mod functions;

use crate::functions::*;

// Determine if the approximation is below an acceptable level of error
fn done(xn: f64, xl: f64, fx: f64) -> bool {
    (xn - xl).abs() + fx.abs() < 1_f64 / 10_f64.powf(8.0)
}

// Bisection Method
fn bisection(f: fn(f64) -> f64, x: (f64, f64), i: i8) -> Option<(f64, i8)> {

    // Find the midpoint
    let c: f64 = (x.1 - x.0) / 2.0 + x.0;

    // Calculate y values
    let fx0: f64 = f(x.0);
    let fc:  f64 = f(c);
    let fx1: f64 = f(x.1);

    // Ensure signs are not the same
    if fx0 * fx1 >= 0.0 {
        println!("The conditions of IVT are not met: f({}) = {} and f({}) = {} are the same sign", x.0, fx0, x.1, fx1);
        return None;
    }

    // Determine new bounds and recurse
    if done(x.0, c, f(c).abs()) {
        Some((c, i))
    } else {
        if i >= std::i8::MAX {
            println!("Function did not converge after {} iterations.", i);
            None
        } else {
            if (fx0 > 0.0) == (fc > 0.0) { bisection(f, (c, x.1), i + 1) } else { bisection(f, (x.0, c), i + 1) }
        }
    }
}

// Newton's Method
fn newton(f: fn(f64) -> f64, ff: fn(f64) -> f64, x: f64, i: i8) -> Option<(f64, i8)> {
    if done(x - f(x) / ff(x), x, f(x)) {
        Some((x - f(x) / ff(x), i))
    } else {
        if i >= std::i8::MAX {
            println!("Function did not coverge after {} iterations.", i);
        } else if ff(x).abs() > 1_f64 / 10_f64.powf(6.0) {
            return newton(f, ff, x - f(x) / ff(x), i + 1);
        } else {
            println!("Unable to approximate: |ff(x)| = |ff({})| < ε = {} < {}", x, ff(x).abs(), 1_f64 / 10_f64.powf(6.0));
        }

        None
    }
}

// Secant Method
fn secant(f: fn(f64) -> f64, x: (f64, f64), i: i8) -> Option<(f64, i8)> {
    if !done(x.1, x.0, f(x.1)) {
        if (f(x.1) - f(x.0)).abs() > 1_f64 / 10_f64.powf(8.0) && i < std::i8::MAX {
            secant(f, (x.1, (x.1 - f(x.1) * ((x.1 - x.0) / (f(x.1) - f(x.0))))), i + 1)
        } else {
            println!("Function did not converge after {} iterations.", i);
            None
        }
    } else {
        Some((x.1 - f(x.1) * ((x.1 - x.0) / (f(x.1) - f(x.0))), i))
    }
}


// Muller's Method
fn muller(f: fn(f64) -> f64, x: (f64, f64, f64), i: i8) -> Option<(f64, i8)> {
    if !done(x.1, x.0, f(x.0)) {
        if i < std::i8::MAX {
            let w: f64 = divdiff!([x.2, x.1]; f) + divdiff!([x.2, x.0]; f) - divdiff!([x.1, x.0]; f);
            let d: f64 = f64::sqrt(w.powf(2.0) - 4.0 * f(x.2) * divdiff!([x.2, x.1, x.0]; f));

            if d.is_nan() {
                println!("x3 is a complex number. Stopping.");
                None
            } else {
                let x3: f64 = x.2 + ((-2.0 * f(x.2)) / (if (w - d).abs() > (w + d).abs() { w - d } else { w + d }));

                muller(f, (x.1, x.2, x3), i + 1)
            }
        } else {
            println!("\nFunction did not converge after {} iterations.", i);
            None
        }
    } else {
        Some((x.2, i))
    }
}

// Inverse quadratic iteration
fn invquad(f: fn(f64) -> f64, x: (f64, f64, f64), i: i8) -> Option<(f64, i8)> {
    if i < std::i8::MAX {
        if !done(x.2, x.1, f(x.1)) {
            let a: f64 = (x.2 - x.1) / (f(x.2) - f(x.1));
            let b: f64 = (1.0 / (f(x.2) - f(x.0))) * (((x.2 - x.1) / (f(x.2) - f(x.1))) - ((x.1 - x.0) / (f(x.1) - f(x.0))));
            let x3: f64 = x.2 - a * f(x.2) + b * f(x.2) * f(x.1);

            if !x3.is_nan() {
                invquad(f, (x.1, x.2, x3), i + 1)
            } else {
                println!("x3 is a complex number. Stopping.");
                None
            }
        } else {
            Some((x.2, i))
        }
    } else {
        println!("Function did not converge after {} iterations.", i);
        None
    }
}

// False Position Method
fn falsepos(f: fn(f64) -> f64, x: (f64, f64), i: i8) -> Option<(f64, i8)> {

    // Find the midpoint
    let c: f64 = (x.0 * f(x.1) - x.1 * f(x.0)) / (f(x.1) - f(x.0));

    // Determine new bounds and recurse
    if done(x.0, c, f(c).abs()) {
        Some((c, i))
    } else {
        if i >= std::i8::MAX {
            println!("Function did not converge after {} iterations.", i);
            None
        } else {
            // Calculate y values
            let fx0: f64 = f(x.0);
            let fc:  f64 = f(c);
            let fx1: f64 = f(x.1);

            // Ensure signs are not the same
            if fx0 * fx1 >= 0.0 {
                println!("The conditions of IVT are not met: f({}) = {} and f({}) = {} are the same sign", x.0, fx0, x.1, fx1);
                return None;
            }

            if (fx0 > 0.0) == (fc > 0.0) { falsepos(f, (c, x.1), i + 1) } else { falsepos(f, (x.0, c), i + 1) }
        }
    }
}

fn main() {

    //-0.247
    println!("\nFunction 1: \n");

    print!("Bisection Method: \t\t");
    
    match bisection(f, (-1.0, 0.0), 0) {
        Some(answer) => {
            println!("{0:.8}", answer.0);
            println!(" - [a, b] = [{}, {}]", -1.0, 0.0);
            println!(" - Iterations: {}", answer.1);
        },
        None => {},
    }

    print!("Newton's Method: \t\t");
    
    match newton(f, ff, -1.0, 0) {
        Some(answer) => {
            println!("{0:.8}", answer.0);
            println!(" - x0 = {}", -1.0);
            println!(" - Iterations: {}", answer.1);
        },
        None => {},
    }

    print!("Secant Method: \t\t\t");
    
    match secant(f, (-1.0, 0.0), 0) {
        Some(answer) => {
            println!("{0:.8}", answer.0);
            println!(" - (x0, x1) = [{}, {}]", -1.0, 0.0);
            println!(" - Iterations: {}", answer.1);
        },
        None => {},
    }

    print!("Muller's Method: \t\t");
    
    match muller(f, (-1.0, -0.5, 0.0), 0) {
        Some(answer) => {
            println!("{0:.8}", answer.0);
            println!(" - (x0, x1, x2) = [{}, {}, {}]", -1.0, -0.5, 0.0);
            println!(" - Iterations: {}", answer.1);
        },
        None => {},
    }

    print!("Inverse Quadratic Iteration: \t");
    
    match invquad(f, (-1.0, -0.5, 0.0), 0) {
        Some(answer) => {
            println!("{0:.8}", answer.0);
            println!(" - (x0, x1, x2) = [{}, {}, {}]", -1.0, -0.5, 0.0);
            println!(" - Iterations: {}", answer.1);
        },
        None => {},
    }

    print!("False Position Method: \t\t");
    
    match falsepos(f, (-1.0, 0.0), 0) {
        Some(answer) => {
            println!("{0:.8}", answer.0);
            println!(" - (x0, x1) = [{}, {}]", -1.0, 0.0);
            println!(" - Iterations: {}", answer.1);
        },
        None => {},
    }

    //2.802
    println!("\nFunction 2: \n");

    print!("Bisection Method: \t\t");
    
    match bisection(f, (2.0, 3.0), 0) {
        Some(answer) => {
            println!("{0:.8}", answer.0);
            println!(" - [a, b] = [{}, {}]", 2.0, 3.0);
            println!(" - Iterations: {}", answer.1);
        },
        None => {},
    }

    print!("Newton's Method: \t\t");
    
    match newton(f, ff, 3.0, 0) {
        Some(answer) => {
            println!("{0:.8}", answer.0);
            println!(" - x0 = {}", 3.0);
            println!(" - Iterations: {}", answer.1);
        },
        None => {},
    }

    print!("Secant Method: \t\t\t");
    
    match secant(f, (3.0, 2.0), 0) {
        Some(answer) => {
            println!("{0:.8}", answer.0);
            println!(" - (x0, x1) = [{}, {}]", 3.0, 2.0);
            println!(" - Iterations: {}", answer.1);
        },
        None => {},
    }

    print!("Muller's Method: \t\t");
    
    match muller(f, (3.0, 2.5, 2.0), 0) {
        Some(answer) => {
            println!("{0:.8}", answer.0);
            println!(" - (x0, x1, x2) = [{}, {}, {}]", 3.0, 2.5, 2.0);
            println!(" - Iterations: {}", answer.1);
        },
        None => {},
    }

    print!("Inverse Quadratic Iteration: \t");
    
    match invquad(f, (3.0, 2.5, 2.0), 0) {
        Some(answer) => {
            println!("{0:.8}", answer.0);
            println!(" - (x0, x1, x2) = [{}, {}, {}]", 3.0, 2.5, 2.0);
            println!(" - Iterations: {}", answer.1);
        },
        None => {},
    }

    print!("False Position Method: \t\t");
    
    match falsepos(f, (2.0, 3.0), 0) {
        Some(answer) => {
            println!("{0:.8}", answer.0);
            println!(" - (x0, x1) = [{}, {}]", 2.0, 3.0);
            println!(" - Iterations: {}", answer.1);
        },
        None => {},
    }

    // -1.0471975542604923248291015625
    println!("\nFunction 3: \n");

    print!("Bisection Method: \t\t");
    
    match bisection(g, (-1.25, -0.75), 90) {
        Some(answer) => {
            println!("{0:.8}", answer.0);
            println!(" - [a, b] = [{}, {}]", 2.0, 3.0);
            println!(" - Iterations: {}", answer.1);
        },
        None => {},
    }

    print!("Newton's Method: \t\t");
    
    match newton(g, gg, -1.75, 0) {
        Some(answer) => {
            println!("{0:.8}", answer.0);
            println!(" - x0 = {}", -1.75);
            println!(" - Iterations: {}", answer.1);
        },
        None => {},
    }

    print!("Secant Method: \t\t\t");
    
    match secant(g, (-1.25, -0.75), 0) {
        Some(answer) => {
            println!("{0:.8}", answer.0);
            println!(" - (x0, x1) = [{}, {}]", -1.25, -0.75);
            println!(" - Iterations: {}", answer.1);
        },
        None => {},
    }

    print!("Muller's Method: \t\t");
    
    match muller(g, (-1.25, -1.0, -0.75), 0) {
        Some(answer) => {
            println!("{0:.8}", answer.0);
            println!(" - (x0, x1, x2) = [{}, {}, {}]", -1.25, -1.0, -0.75);
            println!(" - Iterations: {}", answer.1);
        },
        None => {},
    }

    print!("Inverse Quadratic Iteration: \t");
    
    match invquad(g, (-1.25, -1.0, -1.75), 0) {
        Some(answer) => {
            println!("{0:.8}", answer.0);
            println!(" - (x0, x1, x2) = [{}, {}, {}]", -1.25, -1.0, -0.75);
            println!(" - Iterations: {}", answer.1);
        },
        None => {},
    }

    print!("False Position Method: \t\t");
    
    match falsepos(g, (-1.25, -0.75), 0) {
        Some(answer) => {
            println!("{0:.8}", answer.0);
            println!(" - (x0, x1) = [{}, {}]", -1.25, -0.75);
            println!(" - Iterations: {}", answer.1);
        },
        None => {},
    }
}
