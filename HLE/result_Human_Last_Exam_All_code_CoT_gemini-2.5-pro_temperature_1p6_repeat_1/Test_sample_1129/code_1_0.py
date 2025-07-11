import math

# Step 1 & 2: Define the constants based on the problem and the Lawson criterion.
# The Lawson criterion for a self-sustaining D-T fusion reaction (n*tau) is ~1e20 s/m^3.
# We will use this as an optimistic value for the hypothetical Ti-50 fusion.
lawson_criterion_ntau = 1e20  # units: s * m^-3

# Particle confinement time is given as 'optimistic', 1 second.
confinement_time_tau = 1.0  # units: s

# The reaction chamber is a cube with a side length of 10 cm.
side_length_cm = 10.0  # units: cm

# --- Calculations ---

# Step 3: Calculate the required particle density (n).
# n = (n*tau) / tau
required_density_n = lawson_criterion_ntau / confinement_time_tau

# Step 4: Calculate the volume of the reaction chamber (V).
# First, convert the side length from cm to m.
side_length_m = side_length_cm / 100.0
# Then, calculate the volume of the cube.
volume_m3 = side_length_m ** 3

# Step 5: Calculate the total number of atoms (N).
# N = n * V
total_atoms = required_density_n * volume_m3

# --- Output ---

print("This script calculates the minimum number of atoms for a hypothetical fusion reaction.")
print("-" * 50)
print(f"1. Assumed Lawson Criterion (n*tau): {lawson_criterion_ntau:.1e} s/m^3")
print(f"2. Given Confinement Time (tau): {confinement_time_tau} s")
print(f"3. Given Chamber Side Length: {side_length_cm} cm ({side_length_m} m)")
print("-" * 50)
print(f"First, we calculate the required particle density (n):")
print(f"   n = (n*tau) / tau = {lawson_criterion_ntau:.1e} / {confinement_time_tau} = {required_density_n:.1e} atoms/m^3")
print("\nNext, we calculate the volume of the cubic chamber (V):")
print(f"   V = L^3 = {side_length_m}^3 = {volume_m3:.1e} m^3")
print("\nFinally, we calculate the total number of atoms (N = n * V):")
print(f"   N = {required_density_n:.1e} atoms/m^3 * {volume_m3:.1e} m^3")
print(f"\nThe minimum number of titanium-50 atoms required is: {total_atoms:.1e}")

# The final answer format required by the user.
# The calculation is (1e20 / 1.0) * (0.1**3) = 1e20 * 0.001 = 1e17.
final_answer_val = f"{total_atoms:.1e}"
print(f"\n<<<{final_answer_val}>>>")