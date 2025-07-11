import math

# --- Given Parameters ---
E0 = 8.5  # Initial energy in MeV
R = 8.3   # Total range in cm
x = 4.0   # Distance from the source in cm

# --- Calculation ---
# Calculate the remaining range
R_remaining = R - x

# Calculate the magnitude of the energy loss per centimeter using the derived formula
# |dE/dx| = (2/3) * E0 * R^(-2/3) * (R - x)^(-1/3)
energy_loss = (2/3) * E0 * (R**(-2/3)) * (R_remaining**(-1/3))

# --- Output ---
print("The energy loss per centimetre is calculated based on the empirical Geiger's Law (R ∝ E^(3/2)).")
print("The derived equation for the magnitude of energy loss |dE/dx| at a distance x is:")
print("|dE/dx| = (2/3) * E0 * R^(-2/3) * (R - x)^(-1/3)\n")

print("Substituting the given values into the equation:")
# This line fulfills the requirement to "output each number in the final equation"
print(f"|dE/dx| = (2/3) * {E0} * ({R})**(-2/3) * ({R} - {x})**(-1/3)")
print(f"|dE/dx| = (2/3) * {E0} * ({R})**(-2/3) * ({R_remaining})**(-1/3)\n")

print(f"The calculated energy loss per centimetre at {x} cm is: {energy_loss:.2f} MeV/cm")