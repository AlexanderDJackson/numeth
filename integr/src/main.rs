fn f(x: f64) -> f64 {
    4.0 / (1.0 + x.powf(2.0))
}

fn g(x: f64) -> f64 {
    x.powf(4.0) + x.powf(2.0) + 1.0
}

fn h(x: f64) -> f64 {
    std::f64::consts::E.powf(-x.powf(2.0))
}

fn j(x: f64) -> f64 {
    x.powf(2.0).sin()
}

#[derive(Clone)]
struct Function {
    f: fn(f64) -> f64,
    identifier: &'static str,
    a: f64,
    b: f64,
    n: u16,
    k: u8
}

impl Function {
    pub fn new( f: fn(f64) -> f64, i: &'static str, a: f64, b: f64, n: u16, k: u8) -> Self {
        Self { f: f, identifier: i, a: a, b: b, n: n, k: k }
    }
}

fn prettify(v: Vec<Vec<f64>>) -> String {
    let mut result = String::new();

    for i in v {
        for j in i {
            result.push_str(&(format!("{:.11}\t", j.to_string())));
        }
        result.push('\n');
    }

    String::from(&result[0..result.len() - 2])
}

/*
 * Given: a Function struct
 * Returns: an approximation of the area under f.f from f.a to f.b
 * obtained using a right endpoint riemann sum
 */
fn right(f: &Function) -> f64 {
    let h: f64 = (f.b - f.a) / f.n as f64;
    let mut x = f.a + h;
    let mut sum: f64 = 0.0;

    while x <= f.b {
        sum += (f.f)(x);
        x += h;
    }

    sum * h
}

/*
 * Given: a Function struct
 * Returns: an approximation of the area under f.f from f.a to f.b
 * obtained using a left endpoint riemann sum
 */
fn left(f: &Function) -> f64 {
    let h: f64 = (f.b - f.a) / f.n as f64;
    let mut x = f.a;
    let mut sum: f64 = 0.0;

    while x < f.b {
        sum += (f.f)(x);
        x += h;
    }

    sum * h
}

/*
 * Given: a Function struct
 * Returns: an approximation of the area under f.f from f.a to f.b
 * obtained using a trapezoid riemann sum
 */
fn trapezoid(f: &Function) -> f64 {
    let h: f64 = (f.b - f.a) / f.n as f64;
    let mut x = f.a + h;
    let mut sum: f64 = ((f.f)(f.a) + (f.f)(f.b)) / 2.0;

    for _ in 1..f.n {
        sum += (f.f)(x);
        x += h;
    }

    sum * h
}

/*
 * Given: a Function struct
 * Returns: an approximation of the area under f.f from f.a to f.b
 * obtained using a midpoint endpoint riemann sum
 */
fn midpoint(f: &Function) -> f64 {
    let h: f64 = (f.b - f.a) / f.n as f64;
    let mut x = f.a + h * 0.5;
    let mut sum: f64 = 0.0;

    while x < f.b {
        sum += (f.f)(x);
        x += h;
    }

    sum * h
}

/*
 * Given: a Function struct
 * Returns: an approximation of the area under f.f from f.a to f.b
 * obtained using Simpson's method
 */
fn simpson(f: &Function) -> f64 {
    let h: f64 = (f.b - f.a) / f.n as f64;
    let mut x = f.a + h;
    let mut sum4: f64 = 0.0;

    // Sum up odd values of x to be multiplied by 4
    for _ in 0..(f.n / 2) {
        sum4 += (f.f)(x);
        x += 2.0 * h;
    }

    let mut sum2: f64 = 0.0;
    x = f.a + 2.0 * h;

    // Sum up even values of x to be multiplied by 2
    for _ in 1..(f.n / 2) {
        sum2 += (f.f)(x);
        x += 2.0 * h;
    }

    (h / 3.0) * ((f.f)(f.a) + (f.f)(f.b) + 4.0 * sum4 + 2.0 * sum2)
}

/*
 * Given: a Function struct
 * Returns: A matrix containing successive
 * iterations of richardson's method
 */
fn romberg(f: &Function) -> Vec<Vec<f64>> {
    let mut r = vec![Vec::<f64>::new(); f.k as usize];

    // Generate the approximations with increasing values of n
    for i in 1..(f.k + 1) {
        r[(i - 1) as usize].push(trapezoid(&Function { n: 2_u16.pow(i as u32), ..f.clone() }));
    }

    for i in 1..f.k { // Loop through columns
        for j in i..f.k { // Loop through rows
            
            // Perform richardson's method to improve accuracy
            let new = (4_f64.powf(i as f64) * r[j as usize][(i - 1) as usize] - r[(j - 1) as usize][(i - 1) as usize]) / (4_f64.powf(i as f64) - 1.0);
            r[j as usize].push(new);
        }
    }

    r
}

/*
 * Given: a Function struct
 * Returns: an approximation of the area under the curve
 * using simpson's method with adaptive accuracy
 */
fn adaptive(f: &Function, t: f64) -> f64 {
    let c = (f.a + f.b) / 2.0;

    // Simpson's estimate from a to b
    let sab = ((f.b - f.a) / 6.0) * ((f.f)(f.a) + 4.0 * (f.f)(c) + (f.f)(f.b));

    // Simpson's estimate from a to c
    let sac = ((c - f.a) / 6.0) * ((f.f)(f.a) + 4.0 * (f.f)((f.a + c) / 2.0) + (f.f)(c));

    // Simpson's estimate from c to b
    let scb = ((f.b - c) / 6.0) * ((f.f)(c) + 4.0 * (f.f)((f.b + c) / 2.0) + (f.f)(f.b));

    // If we have the desired accuracy, return the estimate
    if ((sac + scb) - sab).abs() < 15.0 * t {
        sac + scb
    // Otherwise, recurse for the two halves
    } else {
        // The &Function { ..f.clone() } notation describes a function struct that has b = c but
        // otherwise is a copy of f.
        adaptive(&Function { b: c, ..f.clone() }, t) + adaptive(&Function { a: c, ..f.clone() }, t)
    }
}

fn main() {
    /* 
     * Create an array, each element contains:
     * a reference to the Function
     * the identifying string
     * the interval of the estimated integral
     * the number of subsections
     * the the dimensions of the Romberg matrix
     */
    let ident: [Function; 4] = [
        Function::new(f, "4/(1+xˆ2)", 0.0, 1.0, 64, 4),
        Function::new(g, "xˆ4+xˆ2+1", -1.0, 1.0, 64, 4),
        Function::new(h, "exp(-xˆ2)", -2.5, 2.5, 128, 5),
        Function::new(j, "sin(xˆ2)", 0.0, std::f64::consts::PI.sqrt(), 64, 4),
    ];

    let t: f64 = 10_f64.powf(-3.0);

    for i in ident {
        println!("The left endpoint estimate for the Function f(x)={} on the interval [{},{}] with n = {} is {:.11}.",
            i.identifier,
            i.a,
            i.b,
            i.n,
            left(&i));

        println!("The right endpoint estimate for the Function f(x)={} on the interval [{},{}] with n = {} is {:.11}.",
            i.identifier,
            i.a,
            i.b,
            i.n,
            right(&i));

        println!("The trapezoid endpoint estimate for the Function f(x)={} on the interval [{},{}] with n = {} is {:.11}.",
            i.identifier,
            i.a,
            i.b,
            i.n,
            trapezoid(&i));

        println!("The midpoint endpoint estimate for the Function f(x)={} on the interval [{},{}] with n = {} is {:.11}.",
            i.identifier,
            i.a,
            i.b,
            i.n,
            midpoint(&i));

        println!("The Simpson's endpoint estimate for the Function f(x)={} on the interval [{},{}] with n = {} is {:.11}.",
            i.identifier,
            i.a,
            i.b,
            i.n,
            simpson(&i));

        let r = romberg(&i);

        println!("The Romberg algorithm estimate for the Function f(x)={} on the interval [{}, {}] with k = {} is {:.11}.",
            i.identifier,
            i.a,
            i.b,
            i.k,
            r[i.k as usize - 1][i.k as usize - 1]);

        println!("The Romberg matrix is:\n{}", prettify(r));

        println!("The adaptive integration routine for the function f(x)={} on the interval [{}, {}] with tol = {} is {}.\n",
        i.identifier,
        i.a,
        i.b,
        t,
        adaptive(&i, t));
    }
}
