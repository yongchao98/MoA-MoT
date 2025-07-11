import math

# --- Given Information ---
# Initial energy of alpha-particles
E0 = 8.5  # in MeV

# Total range of these particles in standard air
R0 = 8.3  # in cm

# Distance from the source at which to calculate the energy loss
x = 4.0   # in cm

# --- Calculation Plan ---
# 1. The relationship between range (R) and energy (E) is given by Geiger's law: R = a * E^(3/2).
# 2. The energy loss per unit distance, |dE/dx|, can be derived from this law.
#    The derived formula is: |dE/dx| = (2/3) * E₀ * R₀^(-2/3) * (R₀ - x)^(-1/3)
# 3. We will plug the given values into this formula to find the answer.

# --- Step-by-step Calculation ---
# Calculate the remaining range at distance x
remaining_range = R0 - x

# Calculate the energy loss per centimeter using the derived formula
# energy_loss = (2/3) * E₀ * (R₀)^(-2/3) * (R₀ - x)^(-1/3)
energy_loss_per_cm = (2/3) * E0 * math.pow(R0, -2/3) * math.pow(remaining_range, -1/3)

# --- Output the Results ---
print("This script calculates the energy loss per centimeter for an \u03b1-particle in standard air.")
print("\n--- Given Parameters ---")
print(f"Initial Energy (E\u2080): {E0} MeV")
print(f"Total Range (R\u2080): {R0} cm")
print(f"Distance from source (x): {x} cm")

print("\n--- Formula ---")
print("The energy loss |dE/dx| is calculated using the formula derived from Geiger's law:")
print("|dE/dx| = (2/3) * E\u2080 * (R\u2080)^(-2/3) * (R\u2080 - x)^(-1/3)")

print("\n--- Calculation ---")
print(f"Plugging in the values:")
print(f"|dE/dx| = (2/3) * {E0} * ({R0})^(-2/3) * ({R0} - {x})^(-1/3)")
print(f"|dE/dx| = (2/3) * {E0} * ({R0})^(-2/3) * ({remaining_range})^(-1/3)")

print("\n--- Final Answer ---")
print(f"The calculated energy loss per centimetre at {x} cm is: {energy_loss_per_cm:.4f} MeV/cm")

# Output the final numerical answer in the required format
# <<<energy_loss_per_cm>>>