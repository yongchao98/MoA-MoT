import math

#
# Solves for the energy loss per centimeter of an alpha-particle in air.
#

# Step 1: Define initial parameters and determine the constant 'a'
E_initial = 8.5  # Initial energy in MeV
R_total = 8.3  # Total range in cm
distance_x = 4.0  # Distance from the source in cm

# The range-energy relationship is R = a * E^(3/2)
# From this, a = R_total / E_initial^(3/2)
a = R_total / math.pow(E_initial, 1.5)

# Step 2: Calculate the energy at the specified distance
# The remaining range after traveling 'distance_x'
R_remaining = R_total - distance_x

# The energy corresponding to the remaining range
# E = (R / a)^(2/3)
E_at_x = math.pow(R_remaining / a, 2.0 / 3.0)

# Step 3: Calculate the energy loss rate dE/dx
# Differentiating R = a * E^(3/2) gives dR/dE = a * (3/2) * E^(1/2)
# Therefore, the energy loss rate is dE/dx = 1 / (a * (3/2) * E^(1/2))
# dE/dx = (2 / (3 * a)) * E^(-1/2)
energy_loss_per_cm = (2.0 / (3.0 * a)) * math.pow(E_at_x, -0.5)

# Step 4: Display the process and the final answer
print("Calculating the energy loss for an alpha-particle in standard air.")
print(f"Initial Energy (E_initial): {E_initial} MeV")
print(f"Total Range (R_total): {R_total} cm")
print(f"Distance from source (x): {distance_x} cm\n")

print(f"First, we find the energy of the particle after it has traveled {distance_x} cm.")
print(f"The remaining range is R_remaining = {R_total} cm - {distance_x} cm = {R_remaining:.1f} cm.")
print(f"The particle's energy at this point (E_at_x) is {E_at_x:.2f} MeV.\n")

print("The energy loss per centimeter (dE/dx) is then calculated using the formula:")
print("dE/dx = (2 / (3 * a)) * E_at_x^(-1/2)")
print("Plugging in the values:")
print(f"dE/dx = (2 / (3 * {a:.4f})) * {E_at_x:.2f}^(-1/2)")
print(f"dE/dx = {energy_loss_per_cm:.2f} MeV/cm\n")

print(f"The energy loss per centimetre at a distance of {distance_x} cm is {energy_loss_per_cm:.2f} MeV/cm.")

# Final answer formatted as requested
final_answer = round(energy_loss_per_cm, 2)
# <<<0.85>>> is not a valid format as it is not a single character, so just output the number
# It asks for <<<answer content>>> so the number should be inside