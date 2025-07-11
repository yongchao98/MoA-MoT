import math

# We can assume a unit length for the legs of the isosceles right triangle,
# as the final probability is a ratio and is independent of the triangle's actual size.
L = 1.0

# Based on the geometry, we calculate the lengths of the sides.
# The total length of the two legs.
length_legs = 2 * L

# The length of the hypotenuse, calculated using the Pythagorean theorem.
length_hypotenuse = math.sqrt(L**2 + L**2)

# The total perimeter is the sum of the lengths of all sides.
total_perimeter = length_legs + length_hypotenuse

# The escape probability through a side is proportional to its length.
# We calculate the probability of escaping through the hypotenuse (P_hyp).
prob_hypotenuse = length_hypotenuse / total_perimeter

# We calculate the probability of escaping through either of the two legs (P_legs).
prob_legs = length_legs / total_perimeter

# The problem asks for the difference between these two probabilities.
difference = prob_hypotenuse - prob_legs

# The final equation is: Difference = P_hyp - P_legs.
# We will print the values for each component of this equation.

print("The final equation is: Difference = P_hyp - P_legs")
print("The value for each term is calculated below:")

print(f"\n1. Probability of escape through the hypotenuse (P_hyp):")
print(f"P_hyp = {prob_hypotenuse}")

print(f"\n2. Probability of escape through the legs (P_legs):")
print(f"P_legs = {prob_legs}")

print(f"\n3. The final difference is P_hyp - P_legs:")
print(f"Difference = {prob_hypotenuse} - {prob_legs}")
print(f"Difference = {difference}")

# The exact analytical result is 2*sqrt(2) - 3.
exact_difference = 2 * math.sqrt(2) - 3
print(f"\nThe exact value of the difference is 2*sqrt(2) - 3, which is {exact_difference}.")