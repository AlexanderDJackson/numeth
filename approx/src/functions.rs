macro_rules! divdiff {
    () => {};
    ([$left:expr, $right:expr $(,)?]; $f:ident) => {
        ($f($right) - $f($left)) / ($right - $left)
    };
    ([$x:expr $(,)?]; $f:ident) => {
        $f($x)
    };
    ([$first:expr, $second:expr $(, $rest:expr)+ $(,)?]; $f:ident) => {
        divdiff!(@split [$first] [] [$second $(, $rest)+]; $f)
    };
    (@split [$leftmost:expr] [$($middle:expr),*] [$rightmost:expr]; $f:ident) => {
        (divdiff!([$($middle,)* $rightmost]; $f) - divdiff!([$leftmost, $($middle),*]; $f)) / ($rightmost - $leftmost)
    };
    (@split [$leftmost:expr] [$($middle:expr),*] [$next:expr $(, $rest:expr)*]; $f:ident) => {
        divdiff!(@split [$leftmost] [$($middle, )* $next] [$($rest),*]; $f)
    };
}

// Function 1
pub fn f(x: f64) -> f64 {
    x.powf(3.0) - 4.0 * x.powf(2.0) + 3.0 * x + 1.0
}

// Derivative of Function 1
pub fn ff(x: f64) -> f64 {
    3.0 * x.powf(2.0) - 8.0 * x + 3.0
}

// Function 2
pub fn g(x: f64) -> f64 {
    std::f64::consts::E.powf(x) * (1.0 - 2.0 * x.cos())
}

// Derivative of Function 2
pub fn gg(x: f64) -> f64 {
    std::f64::consts::E.powf(x) - 2.0 * std::f64::consts::E.powf(x) * x.cos() + 2.0 * std::f64::consts::E.powf(x) * x.sin()
}
