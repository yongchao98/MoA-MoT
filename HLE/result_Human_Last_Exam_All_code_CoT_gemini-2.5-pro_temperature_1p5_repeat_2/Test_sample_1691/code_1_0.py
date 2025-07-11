import math

# Define the parameters from the integral's denominator polynomial:
# 9.0 * x**5.0 + 5.0 * x**6.0 + 9.0 * x**8.0
# The approximation for small x is dominated by the term with the lowest exponent, p0 * x**n0.
p0 = 9.0
n0 = 5.0

# The analytical approximation for the integral I(epsilon) for small epsilon has the form:
# I(epsilon) ~= C * epsilon**alpha

# The exponent alpha is determined by the dominant power n0.
# alpha = 1/n0 - 1
alpha = 1.0 / n0 - 1.0

# The coefficient C is derived from the parameters and a standard integral result.
# C = p0**(-1/n0) * integral_0^inf(1/(1+u**n0)) du
# The standard integral is equal to (pi/n0) / sin(pi/n0).
integral_val = (math.pi / n0) / math.sin(math.pi / n0)
C = (p0**(-1.0 / n0)) * integral_val

# Output the results, showing how the formula is constructed.
print("The analytical formula that approximates I(epsilon) for small epsilon is derived from its dominant behavior near x=0.")
print("The derived formula structure is: I(epsilon) ~= C * epsilon**alpha")
print("\nBased on the specific integral provided, the detailed formula is:")
print(f"I(epsilon) ~= ( {p0}**(-1.0/{n0}) * (pi/{n0}) / sin(pi/{n0}) ) * epsilon**(1.0/{n0} - 1.0)")

print("\nEvaluating the numerical constants:")
print(f"The exponent alpha = 1.0/{n0} - 1.0 = {alpha}")
print(f"The coefficient C = {p0}**(-1.0/{n0}) * (math.pi/{n0}) / math.sin(math.pi/{n0}) = {C}")

print("\nThus, the final numerical approximation is:")
print(f"I(epsilon) ~= {C:.4f} * epsilon**({alpha})")
