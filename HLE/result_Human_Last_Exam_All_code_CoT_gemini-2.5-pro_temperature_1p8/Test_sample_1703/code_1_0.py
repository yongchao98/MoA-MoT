import math

# The problem asks for the area of triangle T(t) as a function of time t.

# Step 1: Explain the components of the area formula.
# The area A(t) can be expressed as a constant term plus a time-dependent term.
# A(t) = A_midpoint + Delta_A(t)
# A_midpoint is the area when the vertices are at the midpoints of the hexagon sides (t=0).
# Delta_A(t) is the change in area due to the vertices' displacement from the midpoints.

# Step 2: Calculate the coefficients of the formula.
# From geometric derivation, the area as a function of displacement 's' from the midpoint is:
# A(s) = (225 * sqrt(3) / 4) + (3 * sqrt(3) / 4) * s^2
# The constant C1 is the area at s=0. The coefficient C2 is the multiplier for s^2.

R = 10.0
C1_val = (225 * math.sqrt(3)) / 4
C2_val = (3 * math.sqrt(3)) / 4

# Step 3: Define the displacement function.
# The squared displacement of the vertices from the midpoints, s(t)^2, is a periodic
# function of time. It can be expressed compactly as: (5 - abs((t % 10) - 5))^2

# Step 4: Output the final equation with all numbers.
# The instruction requires printing the full equation with numerical values.

print("The area A(t) of the triangle is a periodic function of time.")
print("The formula is composed of a constant base area and a time-dependent part based on the squared displacement of the vertices.")
print("\nThe final equation for the area A(t) is:")
# Note: The 't % 10' in the string represents the modulo operator.
print(f"\nA(t) = {C1_val} + {C2_val} * (5 - abs((t % 10) - 5))^2")

print("\n--- Numerical Values in the Equation ---")
print(f"Base Area (at t=0, 10, 20...): {C1_val}")
print(f"Coefficient for the displacement term: {C2_val}")
print("Periodic displacement function term: (5 - abs((t % 10) - 5))^2")