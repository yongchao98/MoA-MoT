import math

# Based on the step-by-step derivation, the supremum of the angle alpha, M(theta),
# for a small angle theta is found to be approximately M(theta) = (3/4)*theta.

# The problem asks for the limit of M(theta) as theta goes to 0.
# We can represent this relationship in the code.
numerator = 3
denominator = 4

# The final equation is lim_{theta -> 0} M(theta) = L
# From our derivation, M(theta) is proportional to theta, so the limit is 0.
limit_value = 0

print("The analysis of the billiard trajectory for a small angle theta leads to the equation for the supremum M(theta):")
print(f"M(theta) \u2248 ({numerator}/{denominator}) * theta")
print("\nThe numbers in this equation are:")
print(f"Numerator: {numerator}")
print(f"Denominator: {denominator}")
print("\nWe are asked to find the limit of M(theta) as theta approaches 0.")
print(f"Limit = lim_{{theta->0}} ({numerator}/{denominator}) * theta = {limit_value}")
print("\nFinal Answer:")
print(limit_value)