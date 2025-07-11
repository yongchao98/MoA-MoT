import math

# --- Given Parameters ---
# Initial energy of alpha-particles in MeV
E0 = 8.5
# Total range in standard air in cm
R = 8.3
# Distance from the source in cm
x = 4.0

# --- Calculation using the derived formula for energy loss (-dE/dx) ---
# Formula: Energy Loss = (2/3) * (E0 / R) * (1 - x/R)^(-1/3)
term1 = 2 / 3
term2 = E0 / R
term3 = math.pow(1 - (x / R), -1/3)

energy_loss_per_cm = term1 * term2 * term3

# --- Output the final result ---
print(f"To calculate the energy loss per centimetre for α-particles, we use the formula derived from Geiger's rule:")
print(f"Energy Loss = (2/3) * (E₀ / R) * (1 - x / R)⁻¹/³")
print("\nSubstituting the given values:")
print(f"E₀ (Initial Energy) = {E0} MeV")
print(f"R (Total Range) = {R} cm")
print(f"x (Distance) = {x} cm")
print("\nThe equation becomes:")
print(f"Energy Loss = (2/3) * ({E0} MeV / {R} cm) * (1 - {x} cm / {R} cm)⁻¹/³")
print(f"Energy Loss = {term1:.4f} * {term2:.4f} MeV/cm * ({1 - (x / R):.4f})⁻¹/³")
print(f"Energy Loss = {term1:.4f} * {term2:.4f} MeV/cm * {term3:.4f}")
print(f"\nCalculated Energy Loss per Centimetre = {energy_loss_per_cm:.4f} MeV/cm")

# Store the final numerical answer in a variable for easy retrieval.
final_answer = energy_loss_per_cm
# print(f"\n{final_answer}")