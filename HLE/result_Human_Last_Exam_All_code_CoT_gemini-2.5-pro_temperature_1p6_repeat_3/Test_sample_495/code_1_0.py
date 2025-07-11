import math

# --- Given initial values ---
E0 = 8.5  # Initial energy in MeV
R0 = 8.3  # Total range in cm
x = 4.0   # Distance from the source in cm

# --- Step-by-step calculation ---

print("Problem: Calculate the energy loss per centimetre for an α-particle at 4 cm from a source.")
print("-" * 75)
print(f"Initial energy E₀ = {E0} MeV")
print(f"Total range in standard air R₀ = {R0} cm")
print(f"Distance from source x = {x} cm")
print("-" * 75)

# 1. We use the empirical range-energy relationship: R = k * E^(3/2)
# First, we find the constant k using the initial conditions.
k = R0 / (E0 ** 1.5)
print("Step 1: Calculate the proportionality constant 'k' from R = k * E^(3/2)")
print(f"k = R₀ / E₀^(3/2) = {R0} / ({E0}^1.5) = {k:.4f} cm/MeV^(3/2)\n")

# 2. Find the remaining range R' after the particle has traveled x cm.
R_prime = R0 - x
print("Step 2: Calculate the remaining range R' at the specified distance x")
print(f"R' = R₀ - x = {R0} - {x} = {R_prime:.1f} cm\n")

# 3. Calculate the energy E' of the particle at distance x.
# E' corresponds to the energy required to travel the remaining range R'.
# From R' = k * E'^(3/2), we get E' = (R' / k)^(2/3).
E_prime = (R_prime / k) ** (2.0 / 3.0)
print("Step 3: Calculate the energy E' of the particle at that point")
print(f"E' = (R' / k)^(2/3) = ({R_prime:.1f} / {k:.4f})^(2/3) = {E_prime:.3f} MeV\n")

# 4. Calculate the energy loss per centimeter, |dE/dx|.
# The rate of energy loss |dE/dx| is related to the derivative of the range-energy formula.
# dR/dE = (3/2) * k * E^(1/2). The energy loss |dE/dx| = 1 / (dR/dE).
# We evaluate this at E = E'.
energy_loss = 1.0 / ((3.0 / 2.0) * k * (E_prime ** 0.5))

print("Step 4: Calculate the energy loss per centimetre |dE/dx| at E'")
print("From R = k * E^(3/2), we find dR/dE = (3/2) * k * E^(1/2)")
print("The energy loss is |dE/dx| = 1 / (dR/dE)")
print(f"|dE/dx| = 1 / [ (3/2) * {k:.4f} * ({E_prime:.3f}^(1/2)) ]")
print(f"Final calculated energy loss: {energy_loss:.3f} MeV/cm")

<<<0.849>>>