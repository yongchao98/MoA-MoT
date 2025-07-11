import math

# Step 1: Define initial conditions and calculate the constant 'k'
E0 = 8.5  # Initial energy in MeV
R0 = 8.3  # Total range in cm
x = 4.0   # Distance from the source in cm

# Using the relation R = k * E^(3/2)
k = R0 / (E0 ** 1.5)
print(f"Initial Energy (E0): {E0} MeV")
print(f"Total Range (R0): {R0} cm")
print(f"Calculated constant (k): {k:.4f} cm/MeV^(3/2)\n")

# Step 2: Calculate the energy of the particle at x = 4 cm
# The remaining range at x
R_prime = R0 - x
# The energy corresponding to the remaining range
E_at_x = (R_prime / k) ** (2/3)
print(f"Distance from source (x): {x} cm")
print(f"Remaining range at {x} cm (R'): {R_prime:.1f} cm")
print(f"Energy at {x} cm (E(x)): {E_at_x:.4f} MeV\n")

# Step 3 & 4: Calculate the energy loss per cm (dE/dx) at x = 4 cm
# The formula for the magnitude of energy loss is |dE/dx| = 1 / (k * (3/2) * E^(1/2))
energy_loss_per_cm = 1 / (k * 1.5 * math.sqrt(E_at_x))

# Output the final equation with values
print("The energy loss per cm is calculated using the formula: |dE/dx| = 1 / (k * 1.5 * sqrt(E(x)))")
print(f"|dE/dx| = 1 / ({k:.4f} * 1.5 * sqrt({E_at_x:.4f}))")
print(f"|dE/dx| = 1 / ({k:.4f} * 1.5 * {math.sqrt(E_at_x):.4f})")
print(f"|dE/dx| = 1 / ({k * 1.5 * math.sqrt(E_at_x):.4f})")
print(f"\nCalculated energy loss at {x} cm: {energy_loss_per_cm:.4f} MeV/cm")

# Final answer variable for easy extraction
final_answer = energy_loss_per_cm