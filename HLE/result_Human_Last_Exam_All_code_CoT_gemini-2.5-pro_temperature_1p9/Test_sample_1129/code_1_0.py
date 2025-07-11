import math

# --- Step 1: Define Constants and Initial Values ---
# Given values from the problem
temp_celsius = 20.0        # Room temperature in Celsius
chamber_side_cm = 10.0     # Chamber side length in cm
confinement_time_tau = 1.0 # Confinement time in seconds ('optimistic')

# Physical and fusion-specific constants
# Boltzmann constant in eV/K
BOLTZMANN_CONSTANT_EV_K = 8.617333262e-5
# We use the Lawson criterion for D-T fusion as a baseline for the "optimistic" scenario.
# The value is for the triple product (n * T * tau) in units of keV*s/m^3.
LAWSON_CRITERION = 5e21

# --- Step 2: Calculations ---

# Convert temperature from Celsius to Kelvin
temp_kelvin = temp_celsius + 273.15

# Convert temperature from Kelvin to keV
# Energy (eV) = Temperature (K) * Boltzmann Constant (eV/K)
# Energy (keV) = Energy (eV) / 1000
temp_kev = temp_kelvin * BOLTZMANN_CONSTANT_EV_K / 1000

# Convert chamber side length from cm to m
chamber_side_m = chamber_side_cm / 100.0

# Calculate chamber volume in m^3
chamber_volume_m3 = chamber_side_m ** 3

# Calculate the required particle density (n) by rearranging the Lawson criterion:
# n = LAWSON_CRITERION / (T * tau)
particle_density_n = LAWSON_CRITERION / (temp_kev * confinement_time_tau)

# Calculate the total number of atoms (N) required to fill the volume at that density
# N = n * V
total_atoms = particle_density_n * chamber_volume_m3

# --- Step 3: Print the explanation and final result ---
print("To find the minimum number of atoms for sustained fusion, we apply the Lawson criterion.")
print(f"We'll use the optimistic D-T fusion criterion: n * T * τ >= {LAWSON_CRITERION:.0e} keV*s/m³\n")

print(f"1. Convert Temperature to keV:")
print(f"   T = ({temp_celsius}°C + 273.15) * {BOLTZMANN_CONSTANT_EV_K:.6e} eV/K / 1000 = {temp_kev:.6e} keV\n")

print(f"2. Calculate Required Particle Density (n):")
print(f"   n = ({LAWSON_CRITERION:.0e} keV*s/m³) / ({temp_kev:.6e} keV * {confinement_time_tau} s)")
print(f"   n ≈ {particle_density_n:.4e} atoms/m³\n")

print(f"3. Calculate Chamber Volume (V):")
print(f"   V = ({chamber_side_m} m)³ = {chamber_volume_m3} m³\n")

print("4. Calculate Total Required Atoms (N = n * V):")
# The final equation with all the numbers
print(f"   N = {particle_density_n:.4e} atoms/m³ * {chamber_volume_m3} m³")
print(f"   N ≈ {total_atoms:.4e} atoms\n")

# Format the final answer for the required output format
final_answer = f"{total_atoms:.4e}"
print(f"<<<{final_answer}>>>")