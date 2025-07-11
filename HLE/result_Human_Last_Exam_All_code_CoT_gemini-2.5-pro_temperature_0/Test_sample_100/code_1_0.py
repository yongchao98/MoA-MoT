import math

# Based on the step-by-step derivation, the integral evaluates to an expression
# of the form C1 * pi**8 + C2 * pi**2 + C3 * pi + C4.
# The coefficients are derived as follows:
# C1 from the integral of p**7 / (e**p - 1)
# C2 from the sum of integrals of p / (e**p - 1) and p*e**-p / (e**p - 1)
# C3 from the integral of sinh(p/4) / (e**p - 1)
# C4 from the sum of constant terms from the integrals.

# Coefficients from the analytical solution
C1 = 8/15
C2 = 1/3
C3 = -1/2
C4 = 1

# Calculate the terms
term1 = C1 * math.pi**8
term2 = C2 * math.pi**2
term3 = C3 * math.pi
term4 = C4

# The final value of the integral
integral_value = term1 + term2 + term3 + term4

# As requested, we print each number in the final equation.
# The final expression for the integral is: (8/15)*pi**8 + (1/3)*pi**2 - (1/2)*pi + 1
print("The final expression for the integral is of the form: C1 * pi**8 + C2 * pi**2 + C3 * pi + C4")
print("The values of the coefficients are:")
print(f"C1 = {C1}")
print(f"C2 = {C2}")
print(f"C3 = {C3}")
print(f"C4 = {C4}")

print(f"\nFinal Equation: ({C1}) * pi**8 + ({C2}) * pi**2 + ({C3}) * pi + ({C4})")
print(f"\nNumerical value of the integral: {integral_value}")