# About

This project, written by Alex Jackson, implements the following six different root-finding methods:

 - Bisection method
 - Newton’s method
 - Secant method
 - Muller’s method
 - Inverse quadratic iteration
 - False position method (Bracketed secant method)

 Because this was written for a class, I have hard-coded the functions and bounds.

# Running

Ensure you are in the approx project folder.

`cargo run`

# Output

Should look like the following:

```
Function 1: 

Bisection Method: 		-0.24697960
 - [a, b] = [-1, 0]
 - Iterations: 28
Newton's Method: 		-0.24697960
 - x0 = -1
 - Iterations: 5
Secant Method: 			-0.24697960
 - (x0, x1) = [-1, 0]
 - Iterations: 7
Muller's Method: 		-0.24697960
 - (x0, x1, x2) = [-1, -0.5, 0]
 - Iterations: 6
Inverse Quadratic Iteration: 	-0.24697960
 - (x0, x1, x2) = [-1, -0.5, 0]
 - Iterations: 6
False Position Method: 		-0.24697960
 - (x0, x1) = [-1, 0]
 - Iterations: 46

Function 2: 

Bisection Method: 		2.80193773
 - [a, b] = [2, 3]
 - Iterations: 28
Newton's Method: 		2.80193774
 - x0 = 3
 - Iterations: 4
Secant Method: 			2.80193774
 - (x0, x1) = [3, 2]
 - Iterations: 11
Muller's Method: 		1.44504187
 - (x0, x1, x2) = [3, 2.5, 2]
 - Iterations: 7
Inverse Quadratic Iteration: 	2.80193774
 - (x0, x1, x2) = [3, 2.5, 2]
 - Iterations: 19
False Position Method: 		2.80193774
 - (x0, x1) = [2, 3]
 - Iterations: 12

Function 3: 

Bisection Method: 		-1.04719755
 - [a, b] = [2, 3]
 - Iterations: 116
Newton's Method: 		7.33038286
 - x0 = -1.75
 - Iterations: 7
Secant Method: 			-1.04719755
 - (x0, x1) = [-1.25, -0.75]
 - Iterations: 6
Muller's Method: 		-1.04719755
 - (x0, x1, x2) = [-1.25, -1, -0.75]
 - Iterations: 4
Inverse Quadratic Iteration: 	-1.04719755
 - (x0, x1, x2) = [-1.25, -1, -0.75]
 - Iterations: 5
False Position Method: 		-1.04719755
 - (x0, x1) = [-1.25, -0.75]
 - Iterations: 10
```
