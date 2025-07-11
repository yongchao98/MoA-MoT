import math

# Given values
E0 = 8.5  # Initial energy in MeV
R = 8.3   # Total range in cm
x = 4.0   # Distance from the source in cm
n = 1.5   # Geiger's exponent for alpha particles in air

# Calculate the energy loss per centimeter |dE/dx|
# Formula: |dE/dx| = (E0 / (n * R)) * (1 - x/R)**((1-n)/n)
term1 = E0 / (n * R)
base = 1 - (x / R)
exponent = (1 - n) / n
term2 = math.pow(base, exponent)
energy_loss_per_cm = term1 * term2

# Print the equation with the values
print("Calculation of energy loss per centimeter |dE/dx|:")
print(f"|dE/dx| = ({E0} MeV / ({n} * {R} cm)) * (1 - {x} cm / {R} cm) ** (({1} - {n}) / {n})")
print(f"|dE/dx| = {energy_loss_per_cm:.3f} MeV/cm")
