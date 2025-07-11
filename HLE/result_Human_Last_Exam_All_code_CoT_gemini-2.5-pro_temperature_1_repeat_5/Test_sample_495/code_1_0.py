import math

# Initial values provided in the problem
E0 = 8.5  # Initial energy in MeV
R0 = 8.3  # Total range in standard air in cm
x = 4.0   # Distance from the source in cm

print("This script calculates the energy loss per centimetre for an alpha-particle in air.")
print("-" * 60)

# Step 1: Calculate the residual range of the alpha-particle
R_res = R0 - x
print("Step 1: Calculate the residual range (R_res)")
print(f"The particle has traveled x = {x:.1f} cm from a total possible range of R0 = {R0:.1f} cm.")
print(f"The residual range is R_res = R0 - x = {R0:.1f} - {x:.1f} = {R_res:.1f} cm.")
print("-" * 60)

# Step 2: Calculate the residual energy (E_x) at the distance x
# This is derived from the empirical relation R = a * E^(3/2)
# E_x = E0 * (R_res / R0)^(2/3)
E_x = E0 * math.pow(R_res / R0, 2/3)
print("Step 2: Calculate the residual energy (E_x) at this point")
print("Using the formula E_x = E0 * (R_res / R0)^(2/3):")
print(f"E_x = {E0:.1f} * ({R_res:.1f} / {R0:.1f})^(2/3)")
print(f"E_x = {E_x:.3f} MeV")
print("-" * 60)

# Step 3: Calculate the energy loss per centimeter (|dE/dx|)
# This is derived from dE/dx = 1 / (dR/dE)
# It simplifies to |dE/dx| = (2/3) * (E_x / R_res)
dEdx = (2/3) * (E_x / R_res)
print("Step 3: Calculate the energy loss per centimeter (|dE/dx|)")
print("Using the formula |dE/dx| = (2/3) * E_x / R_res:")
print(f"|dE/dx| = (2/3) * {E_x:.3f} MeV / {R_res:.1f} cm")
print(f"|dE/dx| = {dEdx:.3f} MeV/cm")
print("-" * 60)

print(f"Final Answer: The energy loss per centimetre for the α-particles at a distance of {x:.1f} cm is {dEdx:.3f} MeV/cm.")