import math

# Step 1: Explain the problem and the derived formula.
print("This problem is about finding the shape of a given volume that maximizes the gravitational field at a point.")
print("The optimal shape is bounded by the surface r = R_max * sqrt(cos(θ)) in spherical coordinates.")
print("The furthest point on this surface from the point of interest (A) is R_max.")
print("\nThe volume (V) of this shape is given by the formula: V = (4 * pi / 15) * R_max^3")

# Step 2: State the given values and the equation to be solved.
volume = 1.0
print(f"\nWe are given V = {volume} cubic meter.")
print("To find R_max, we solve the equation: 1 = (4 * pi / 15) * R_max^3")
print("Rearranging for R_max gives: R_max = (15 / (4 * pi))^(1/3)")

# Step 3: Define the constants and perform the calculation.
print("\nPerforming the calculation:")
numerator = 15.0
denominator_factor = 4.0
pi = math.pi
power = 1.0 / 3.0

# Calculate R_max
r_max = (numerator / (denominator_factor * pi))**power

# Step 4: Display the final equation with all numbers and the result.
print(f"The final equation is: R_max = ({numerator} / ({denominator_factor} * {pi}))^({power})")
print(f"The furthest point on the surface of the playdough is {r_max} meters from point A.")
