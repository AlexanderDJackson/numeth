fn f(t: f64, y: f64) -> f64 {
    t.powf(2.0) * y.sin() + y * (std::f64::consts::PI * t).sin()
}
//h = 0.1: 3.2914688
//h = 0.001: 3.2914688
fn euler(t: f64, t_final: f64, mut y: f64, h: f64) -> f64 {
    let mut i: f64 = 0.0;
    while i < ((t_final - t) / h) {
        y = y + h * f(t + i * h, y);
        i += 1.0;
    }

    y
}

#[allow(dead_code)]
fn pc(mut t: f64, t_final: f64, mut y: f64, h: f64, c: f64) -> f64 {
    while t <= t_final {
        y = y + h * (c * f(t, y) + (1.0 - c) * f(t + (1.0 / (2.0 * (1.0 - c))) * h, y + (1.0 / (2.0 * (1.0 - c))) * h * f(t, y)));
        t += h;
    }

    y
}

fn rkf4(t: f64, y: f64, h: f64) -> f64 {
    let mut k: [f64; 6] = [0.0; 6];

    k[0] = f(t, y);
    k[1] = f(t + (h / 5.0), y + h * k[0] / 5.0);
    k[2] = f(t + 3.0 * h / 10.0, y + 3.0 * h * k[1] / 40.0 + 9.0 * h * k[1] / 40.0);
    k[3] = f(t + 4.0 * h / 5.0, y + 44.0 * h * k[0] / 45.0 - 56.0 * h * k[1] / 15.0 + 32.0 * h * k[2] / 9.0);
    k[4] = f(t + 8.0 * h / 9.0, y + 19372.0 * h * k[0] / 6561.0 - 25360.0 * h * k[1] / 2187.0 + 64448.0 * h * k[2] / 6561.0 - 212.0 * h * k[3] / 729.0);
    k[5] = f(t + h, y + 9017.0 * h * k[0] / 3168.0 - 355.0 * h * k[1] / 33.0 + 46732.0 * h * k[3] / 5247.0 + 49.0 * h * k[3] / 176.0 - 5103.0 * h * k[4] / 18656.0);

    y + h * (35.0 * k[0] / 384.0 + 500.0 * k[2] / 1113.0 + 125.0 * k[3] / 192.0 - 2187.0 * k[4] / 6784.0 + 11.0 * k[5] / 84.0)
}

fn rkf5(t: f64, y: f64, h: f64) -> f64 {
    let mut k: [f64; 6] = [0.0; 6];

    k[0] = f(t, y);
    k[1] = f(t + h / 5.0, y + h * k[0] / 5.0);
    k[2] = f(t + 3.0 * h / 10.0, y + 3.0 * h * k[0] / 40.0 + 9.0 * h * k[1] / 40.0);
    k[3] = f(t + 4.0 * h / 5.0, y + 44.0 * h * k[0] / 45.0 - 56.0 * h * k[1] / 15.0 + 32.0 * h * k[2] / 9.0);
    k[4] = f(t + 8.0 * h / 9.0, y + 19372.0 * h * k[0] / 6561.0 - 25360.0 * h * k[1] / 2187.0 + 64448.0 * h * k[2] / 6561.0 - 212.0 * h * k[3] / 729.0);
    k[5] = f(t + h, y + 9017.0 * h * k[0] / 3168.0 - 355.0 * h * k[1] / 33.0 + 46732.0 * h * k[2] / 5247.0 + 49.0 * h * k[3] / 176.0 - 5103.0 * h * k[4] / 18656.0);

    y + h * (35.0 * k[0] / 384.0 + 500.0 * k[2] / 1113.0 + 125.0 * k[3] / 192.0 - 2187.0 * k[4] / 6784.0 + 11.0 * k[5] / 84.0)
}

fn rkf45(t: f64, t_final: f64, y: f64, h: f64, e: f64) -> f64 {
    if t + h > t_final || t >= t_final {
        y
    } else {
        let y1 = rkf4(t, y, h);
        let y2 = rkf5(t, y, h);

        let error = (y1 - y2).abs();

        if 0.25 * h * e > error {
            rkf45(t + h, t_final, y2, 2.0 * h, e)
        } else if error > h * e {
            rkf45(t, t_final, y, h / 2.0, e)
        } else {
            rkf45(t + h, t_final, y2, h, e)
        }
    }
}

fn main() {
    println!("Method 1: Euler's Method\n------------------------");
    println!("f({}, {}) = {:.8} with h = {}", 1, 3, euler(1.0, 3.0, 1.0, 0.1), 0.1);
    println!("f({}, {}) = {:.8} with h = {}", 1, 3, euler(1.0, 3.0, 1.0, 0.001), 0.001); //3.2914688
    println!("\nMethod 2: RFK45\n---------------");
    println!("f({}, {}) = {:.8} with h = {}", 1, 3, rkf45(1.0, 3.0, 1.0, 0.1, 10_f64.powf(-6.0)), 0.1);
    println!("f({}, {}) = {:.8} with h = {}", 1, 3, rkf45(1.0, 3.0, 1.0, 0.001, 10_f64.powf(-6.0)), 0.001);
}