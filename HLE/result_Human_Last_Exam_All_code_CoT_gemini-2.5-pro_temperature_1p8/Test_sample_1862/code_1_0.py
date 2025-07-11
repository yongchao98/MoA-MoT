import math

# --- Given Parameters ---

# Case 1: Full spherical particle
diameter1_nm = 50
time1_s = 20

# Case 2: Quarter-sphere particle
diameter2_nm = 250

# --- Calculations ---

# Step 1: Calculate the radii from the diameters
radius1_nm = diameter1_nm / 2
radius2_nm = diameter2_nm / 2

# Step 2: Apply the diffusion time relationship
# The lithiation time (t) is proportional to the square of the radius (r): t = k * r^2
# Therefore, t2 = t1 * (r2 / r1)^2
time2_s = time1_s * (radius2_nm / radius1_nm)**2

# --- Output the results ---

print("This problem is solved based on the principle of diffusion.")
print("The time required for lithiation (t) is proportional to the square of the particle's radius (r).\n")
print(f"The relationship is: t2 = t1 * (r2 / r1)^2\n")

print("--- Given Values ---")
print(f"Initial Radius (r1): {diameter1_nm} nm / 2 = {radius1_nm:.1f} nm")
print(f"Initial Time (t1): {time1_s} s")
print(f"New Radius (r2): {diameter2_nm} nm / 2 = {radius2_nm:.1f} nm\n")


print("--- Final Equation and Calculation ---")
# The final equation with each number explicitly shown
print(f"Lithiation Time (t2) = {time1_s} s * ({radius2_nm:.1f} nm / {radius1_nm:.1f} nm)^2")
print(f"Lithiation Time (t2) = {time1_s} s * ({radius2_nm / radius1_nm:.1f})^2")
print(f"Lithiation Time (t2) = {time1_s} s * {math.pow(radius2_nm / radius1_nm, 2):.1f}")
print(f"Lithiation Time (t2) = {time2_s:.0f} s\n")

print(f"The calculated lithiation time for the 250 nm quarter-sphere particle is {time2_s:.0f} seconds.")