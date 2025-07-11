import math

# Define the constants from the problem
R = 10.0  # Radius of the circle
s = R     # Side length of the inscribed regular hexagon
v = 1.0   # Speed of the vertices

# The derived formula for the area is of the form: A(t) = C * (t^2 + K)
# Let's calculate the constants C and K based on the derivation.

# C = (3 * sqrt(3)) / 4
C_numerator = 3
C_denominator = 4
C_value = (C_numerator * math.sqrt(3)) / C_denominator

# K = (3/4) * s^2
K_numerator = 3
K_denominator = 4
K_value = (K_numerator / K_denominator) * s**2

# Print the final result, showing each number in the equation.
print("The area of the triangle T(t) as a function of time t is:")
print(f"Area(t) = ({C_numerator} * sqrt(3) / {C_denominator}) * (t^2 + {int(K_value)})")
print("\nWhere:")
print(f"The coefficient outside the parenthesis is (3 * sqrt(3) / 4) ≈ {C_value:.4f}")
print(f"The constant added to t^2 is (3/4) * s^2 = {K_value}")
