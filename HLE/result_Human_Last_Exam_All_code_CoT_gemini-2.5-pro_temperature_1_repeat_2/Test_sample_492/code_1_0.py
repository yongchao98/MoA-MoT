import math

# This script calculates the exact time 'c' of the emergence of the giant
# connected component in the specified dynamic random graph model.

# The condition for the emergence is given by the Molloy-Reed criterion:
# E[D(D-1)] / E[D] = 1
# where D is the degree of a randomly chosen vertex.

# Based on the model's properties, we derived the following moments as a function of time t:
# E[D] = t^2 / 3
# E[D(D-1)] = 2 * t^4 / 15

print("The critical time 'c' is found by solving the equation:")
print("(2 * c^4 / 15) / (c^2 / 3) = 1")

print("\nSimplifying the equation leads to:")
print("2 * c^2 / 5 = 1")
print("c^2 = 5 / 2")

# Define the numbers in the final equation c^2 = numerator / denominator
numerator = 5
denominator = 2

# Calculate c
c_squared = numerator / denominator
c = math.sqrt(c_squared)

print(f"\nThe final equation is c^2 = {numerator}/{denominator}.")
print(f"The numbers in this equation are numerator={numerator} and denominator={denominator}.")
print(f"\nThe exact value of c is the square root of ({numerator}/{denominator}).")
print(f"c = sqrt({c_squared})")
print(f"The numerical value of c is approximately: {c}")