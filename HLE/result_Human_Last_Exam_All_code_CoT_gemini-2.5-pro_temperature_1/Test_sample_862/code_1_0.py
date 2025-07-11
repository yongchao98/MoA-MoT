# The problem reduces to finding the minimum of the first non-zero eigenvalue
# of a Sturm-Liouville problem. Through variational methods, we find that
# the wave number k of the extremal function satisfies two equations.

# Let X = k * delta and Y = k * (pi - delta). The equations are:
# tan(Y) = -1/3 * tan(X)
# cos(X)^2 = 1/3 * cos(Y)^2
# Solving this system leads to a minimal wave number k.

# The solution for k gives the following numerator and denominator.
k_numerator = 5
k_denominator = 6

# The smallest possible eigenvalue lambda_min is k^2.
lambda_min_numerator = k_numerator**2
lambda_min_denominator = k_denominator**2

# The constant C is the inverse of this minimal eigenvalue.
C_numerator = lambda_min_denominator
C_denominator = lambda_min_numerator

# Calculate the final value of C.
C = C_numerator / C_denominator

print(f"The minimal value for the wave number k is {k_numerator}/{k_denominator}.")
print(f"This gives a minimal eigenvalue of lambda_min = k^2 = ({k_numerator}/{k_denominator})^2 = {lambda_min_numerator}/{lambda_min_denominator}.")
print(f"The constant C is the inverse of the minimal eigenvalue: C = 1 / lambda_min.")
print(f"Therefore, C = ({k_denominator}/{k_numerator})^2 = {C_numerator}/{C_denominator}.")
print(f"The numerical value of C is {C}.")
