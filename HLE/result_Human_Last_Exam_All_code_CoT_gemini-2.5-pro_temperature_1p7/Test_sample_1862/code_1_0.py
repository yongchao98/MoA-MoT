# Plan:
# 1. Define the known variables from the problem statement for both cases.
# 2. Calculate the radius for each particle from its diameter.
# 3. Apply the diffusion scaling law (t ∝ r²) to calculate the new lithiation time.
# 4. Print the reasoning, the formula with the values substituted, and the final result.

# Case 1: Full sphere
d1_nm = 50.0  # Diameter in nm
t1_s = 20.0   # Lithiation time in seconds

# Case 2: Quarter sphere
d2_nm = 250.0 # Diameter in nm

# Step 1: Calculate the radii
r1_nm = d1_nm / 2.0
r2_nm = d2_nm / 2.0

# Step 2: Apply the diffusion time scaling law t₂ = t₁ * (r₂ / r₁)²
t2_s = t1_s * (r2_nm / r1_nm)**2

# Step 3: Print the explanation and the result
print("The lithiation time is limited by solid-state diffusion, where time (t) is proportional to the square of the particle radius (r).")
print("The relationship is: t₂ = t₁ * (r₂ / r₁)²\n")
print("Given values:")
print(f"t₁ (initial time) = {t1_s} s")
print(f"r₁ (initial radius) = {d1_nm} nm / 2 = {r1_nm} nm")
print(f"r₂ (new radius) = {d2_nm} nm / 2 = {r2_nm} nm\n")
print("Calculating the new lithiation time (t₂):")
print(f"t₂ = {t1_s} s * ({r2_nm} nm / {r1_nm} nm)²")
print(f"t₂ = {t1_s} s * ({r2_nm / r1_nm})²")
print(f"t₂ = {t1_s} s * ({(r2_nm / r1_nm)**2})")
print(f"t₂ = {t2_s} s\n")

print(f"The calculated lithiation time for the 250 nm quarter-sphere particle is {t2_s} seconds.")
print("<<<500.0>>>")