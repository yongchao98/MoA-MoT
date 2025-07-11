import math

# --- Given parameters for the first case (full sphere) ---
t1 = 20  # Lithiation time in seconds
d1 = 50  # Diameter in nm

# --- Given parameters for the second case (quarter sphere) ---
d2 = 250  # Diameter in nm

# --- Step 1: Calculate the radii ---
r1 = d1 / 2
r2 = d2 / 2

# --- Step 2: Apply the diffusion time scaling law ---
# The time 't' for complete lithiation is governed by diffusion and scales with the square
# of the characteristic length (the radius 'r'). So, t ∝ r².
# The relationship is t2 / t1 = (r2 / r1)²
# Therefore, t2 = t1 * (r2 / r1)²
t2 = t1 * (r2 / r1)**2

# --- Step 3: Print the explanation and the final calculation ---
print("The lithiation process is limited by solid-state diffusion.")
print("The time required for complete lithiation (t) is proportional to the square of the particle's radius (r).")
print("The scaling law is: t₂ = t₁ * (r₂ / r₁)²")
print("\nGiven values:")
print(f"Initial time (t₁): {t1} s")
print(f"Initial radius (r₁): {d1}/2 = {r1} nm")
print(f"New radius (r₂): {d2}/2 = {r2} nm")
print("\nCalculation:")
# To show the calculation with intermediate steps as requested
ratio_squared = (r2 / r1)**2
print(f"t₂ = {t1} s * ({r2} nm / {r1} nm)²")
print(f"t₂ = {t1} s * ({int(r2/r1)})²")
print(f"t₂ = {t1} s * {int(ratio_squared)}")
print(f"t₂ = {int(t2)} s")

print(f"\nThe calculated lithiation time for the 250 nm quarter-sphere particle is {int(t2)} seconds.")

# Final answer in the specified format
# print(f'<<<{int(t2)}>>>')