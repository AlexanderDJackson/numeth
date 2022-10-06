#![feature(trace_macros)]

trace_macros!(true);

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

fn main() {

    let x = 4.0;
    
    macro_rules! lagrange {
        // Nothing
        () => {};

        // One point
        ((($x:literal, $y:literal))) => { $y };
        ($x:literal, $y:literal) => { lagrange!((($x, $y))) };
        (($x:literal, $y:literal)) => { lagrange!((($x, $y))) };

        // First term
        ((($x:literal, $y:literal), $(($i:literal, $j:literal)$(,)? )+)) => {
            $y $(* (x - $i) / ($x - $i))+ + lagrange!(($(($i, $j)),+), (($x, $y)))
        };

        // Nth term
        ((($x:literal, $y:literal)$(,)? $(($i:literal, $j:literal)$(,)? )+), (($a:literal, $b:literal), $(($c:literal, $d:literal)),+)) => {
            $y $(* (x - $i) / ($x - $i))+ * (x - $a) / ($x - $a) $(* (x - $c) / ($x - $d))* + lagrange!(($(($i, $j)),+), (($x, $y), ($a, $b), $(($c, $d)),+))
        };

        // Second term
        ((($x:literal, $y:literal), $(($i:literal, $j:literal)$(,)? )+), (($a:literal, $b:literal))) => {
            $y $(* (x - $i) / ($x - $i))+ * (x - $a) / ($x - $a) + lagrange!(($(($i, $j)),+), (($x, $y), ($a, $b)))
        };

        // Last term
        ((($x:literal, $y:literal)), (($a:literal, $b:literal), $(($c:literal, $d:literal)),+)) => {
            $y * (x - $a) / ($x - $a) $(* (x - $c) / ($x - $c))*
        };
    }
    
    macro_rules! l {
        // Nothing
        () => {};

        // One point
        ((($x:literal, $y:literal))) => { $y };
        ($x:literal, $y:literal) => { lagrange!((($x, $y))) };
        (($x:literal, $y:literal)) => { lagrange!((($x, $y))) };

        // First term
        ((($x:literal, $y:literal), $(($i:literal, $j:literal)$(,)? )+)) => {
            $y $(* (x - $i) / ($x - $i))+ + lagrange!(($(($i, $j)),+), (($x, $y)))
        };

        // Nth term
        ((($x:literal, $y:literal)$(,)? $(($i:literal, $j:literal)$(,)? )+), (($a:literal, $b:literal), $(($c:literal, $d:literal)),+)) => {
            $y $(* (x - $i) / ($x - $i))+ * (x - $a) / ($x - $a) $(* (x - $c) / ($x - $d))* + lagrange!(($(($i, $j)),+), (($x, $y), ($a, $b), $(($c, $d)),+))
        };

        // Last term
        ((($x:literal, $y:literal)), (($a:literal, $b:literal), $(($c:literal, $d:literal)),+)) => {
            $y * (x - $a) / ($x - $a) $(* (x - $c) / ($x - $c))*
        };
    }

    //assert!(1.0 * (x - 2.0) / (1.0 - 2.0) * (x - 3.0) / (1.0 - 3.0) + 4.0 * (x - 1.0) / (2.0 - 1.0) * (x - 3.0) / (2.0 - 3.0) + 9.0 * (x - 1.0) / (3.0 - 1.0) * (x - 2.0) / (3.0 - 2.0), lagrange!(((1.0, 1.0), (2.0, 4.0), (3.0, 9.0))));
    println!("{}", l!(((1.0, 1.0), (2.0, 8.0), (3.0, 27.0), (4.0, 64.0))));
}
