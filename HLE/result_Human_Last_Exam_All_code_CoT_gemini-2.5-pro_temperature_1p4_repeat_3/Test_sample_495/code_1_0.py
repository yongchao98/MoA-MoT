import math

# --- Given values ---
# Initial energy of α-particles in MeV
E0 = 8.5
# Total range in standard air in cm
R0 = 8.3
# Distance from the source in cm
x = 4.0

# --- Calculation ---
# 1. Calculate the remaining range of the α-particles
R_rem = R0 - x

# 2. Calculate the energy loss per centimeter using the formula derived from Geiger's rule.
# The formula is: |dE/dx| = (2/3) * E₀ * R₀^(-2/3) * R_rem^(-1/3)
energy_loss = (2/3) * E0 * (R0**(-2/3)) * (R_rem**(-1/3))

# --- Output ---
print("Calculation of energy loss for α-particles in standard air.")
print("-" * 50)
print(f"Initial Energy (E₀): {E0} MeV")
print(f"Total Range (R₀): {R0} cm")
print(f"Distance from source (x): {x} cm")
print(f"Remaining Range (R_rem = R₀ - x): {R0} - {x} = {R_rem} cm")
print("-" * 50)
print("The final equation for the energy loss per centimeter is:")
# Print the equation with all the numbers as requested
print(f"|dE/dx| = (2/3) * {E0} * ({R0})**(-2/3) * ({R_rem})**(-1/3)")
print(f"\nResult:")
print(f"The energy loss per centimetre at {x} cm is {energy_loss:.4f} MeV/cm.")
