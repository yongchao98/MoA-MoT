import math

# The problem asks for the smallest possible value k such that for any smooth
# Riemannian metric on S^2 with total area 4*pi, the second nonzero eigenvalue
# lambda_2 is always < k. This value k is the supremum of lambda_2 over all
# such metrics.

# According to a theorem by Romain Petrides (2014), for any smooth Riemannian
# metric g on the sphere S^2, the second nonzero eigenvalue lambda_2 satisfies
# the inequality:
# lambda_2(g) * Area(g) <= 8 * pi.
# The equality is achieved for the standard round sphere.

# The given area is 4*pi.
area = 4 * math.pi

# The constant factor in the numerator of the inequality is 8*pi.
bound_numerator = 8 * math.pi

# We can calculate the supremum of lambda_2 by rearranging the inequality:
# k = sup(lambda_2) = (8 * pi) / Area
k = bound_numerator / area

# We print the steps of the calculation.
# The equation is k = 8*pi / 4*pi
numerator_val_in_pi = 8
denominator_val_in_pi = 4

print("Based on Petrides' theorem, the supremum of the second nonzero eigenvalue λ₂ on a sphere is bounded:")
print("λ₂ * Area ≤ 8π")
print("\nGiven Area = 4π, we can find the value of k:")
print(f"k = ({numerator_val_in_pi} * π) / ({denominator_val_in_pi} * π)")
print(f"k = {k}")
