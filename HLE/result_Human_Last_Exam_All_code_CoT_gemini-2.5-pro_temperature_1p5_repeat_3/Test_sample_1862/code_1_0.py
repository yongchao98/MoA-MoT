import math

# Initial conditions for the first particle (full sphere)
d1_nm = 50  # diameter in nm
t1_s = 20   # lithiation time in seconds

# Conditions for the second particle (quarter sphere)
d2_nm = 250 # diameter in nm

# Step 1: Calculate the radii from the diameters
r1_nm = d1_nm / 2
r2_nm = d2_nm / 2

# Step 2: Calculate the new lithiation time using the diffusion time relationship
# t2 = t1 * (r2 / r1)^2
t2_s = t1_s * (r2_nm / r1_nm)**2

# Step 3: Print the explanation and the calculation step-by-step
print("The lithiation process is limited by solid-state diffusion.")
print("The time for diffusion (t) is proportional to the square of the diffusion length (radius, r): t ∝ r².")
print("Therefore, the new lithiation time (t2) can be calculated using the formula: t2 = t1 * (r2 / r1)²\n")

print("Given values:")
print(f"t1 = {t1_s} s")
print(f"d1 = {d1_nm} nm => r1 = {r1_nm} nm")
print(f"d2 = {d2_nm} nm => r2 = {r2_nm} nm\n")

print("Calculation:")
# Using an f-string to output the equation with the numbers
print(f"t2 = {t1_s} s * ({r2_nm} nm / {r1_nm} nm)²")
print(f"t2 = {t1_s} s * ({r2_nm / r1_nm})²")
print(f"t2 = {t1_s} s * {(r2_nm / r1_nm)**2}")
print(f"t2 = {t2_s} s\n")

print(f"The final calculated lithiation time for the quarter-sphere particle is {t2_s} seconds.")
<<<500.0>>>