import math

# --- Given parameters for the first case ---
d1_nm = 50.0  # Diameter of the first particle in nm
t1_s = 20.0   # Lithiation time for the first particle in seconds

# --- Given parameters for the second case ---
d2_nm = 250.0 # Diameter of the second particle in nm

# --- Step 1: Calculate the radii ---
r1_nm = d1_nm / 2
r2_nm = d2_nm / 2

print("Step 1: Calculate radii for both particles.")
print(f"Radius of the first particle (r1) = {d1_nm} nm / 2 = {r1_nm} nm")
print(f"Radius of the second particle (r2) = {d2_nm} nm / 2 = {r2_nm} nm\n")

# --- Step 2: Apply the diffusion time relationship ---
# The lithiation time (t) is proportional to the square of the radius (r), so t_new = t_old * (r_new / r_old)^2.
# The shape being a "quarter sphere" does not change the maximum diffusion length, which is its radius.
ratio_r_squared = (r2_nm / r1_nm)**2

# --- Step 3: Calculate the new lithiation time ---
t2_s = t1_s * ratio_r_squared

print("Step 2: Calculate the new lithiation time (t2).")
print("The formula is: t2 = t1 * (r2 / r1)^2")
print(f"t2 = {t1_s} s * ({r2_nm} nm / {r1_nm} nm)^2")
print(f"t2 = {t1_s} s * ({r2_nm / r1_nm})^2")
print(f"t2 = {t1_s} s * {ratio_r_squared}")
print(f"Final calculated lithiation time (t2) = {t2_s} seconds")

# The final answer in the required format will be printed by the user, but this represents the numerical result.
# For the purpose of this script, we print the final value.
# print(f"\n<<< {t2_s} >>>")