from fractions import Fraction

# The invariant density has the form rho(x) = a + b * x^(1/2)
# The unnormalized form is proportional to 1 + 2*sqrt(x).
# We calculate the normalization constant C such that the integral of C*(1+2*sqrt(x)) from 0 to 1 is 1.
# C * integral(1+2*x^(1/2))dx = C * [x + 4/3*x^(3/2)]_0^1 = C * (1 + 4/3) = C * 7/3
# C * 7/3 = 1  => C = 3/7
C = Fraction(3, 7)
a_unnormalized = 1
b_unnormalized = 2

# The coefficients of the normalized density
a = C * a_unnormalized
b = C * b_unnormalized

print("The normalised density of the invariant measure is of the form: rho(x) = a + b * x^(1/2)")
print("The coefficients are:")
print(f"a = {a}")
print(f"b = {b}")
print("\nThe final equation is:")
print(f"rho(x) = {a} + {b} * x^(1/2)")
