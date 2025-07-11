import math

# Given values
E_initial = 8.5  # MeV
R_total = 8.3    # cm
x = 4.0        # cm

# Formula for energy loss per cm: -dE/dx = (2 * E_initial) / (3 * R_total) * (1 - x / R_total)^(-1/3)
term1 = (2 * E_initial) / (3 * R_total)
term2 = (1 - x / R_total)
exponent = -1.0 / 3.0

energy_loss = term1 * (term2 ** exponent)

print("Calculating the energy loss per centimeter (-dE/dx) for an alpha-particle.")
print(f"Initial Energy (E_initial): {E_initial} MeV")
print(f"Total Range (R_total): {R_total} cm")
print(f"Distance from source (x): {x} cm")
print("\nUsing the formula: -dE/dx = (2 * E_initial) / (3 * R_total) * (1 - x/R_total)^(-1/3)")
print("\nCalculation steps:")
print(f"-dE/dx = (2 * {E_initial}) / (3 * {R_total}) * (1 - {x}/{R_total})^(-1/3)")
print(f"-dE/dx = ({2 * E_initial}) / ({3 * R_total}) * ({1 - x/R_total})^(-1/3)")
print(f"-dE/dx = {term1} * {term2**exponent}")
print(f"-dE/dx = {energy_loss:.4f} MeV/cm")

print(f"\nThus, the energy loss per centimetre for these α-particles at a distance of {x} cm from the source is approximately {energy_loss:.2f} MeV/cm.")
