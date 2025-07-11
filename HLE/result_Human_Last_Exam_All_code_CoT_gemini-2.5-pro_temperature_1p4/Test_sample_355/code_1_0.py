import math

# Step 0: Define constants and assumptions for the calculation.

# Telescope diameter in meters
D_m = 1.0
# Star's apparent B-band magnitude (assuming M_B in problem meant m_B)
m_B = 20.0
# Desired magnitude error
delta_m = 0.01
# Zero-point flux density for B-band (a standard value in photons/s/cm^2/Angstrom)
f0_lambda = 1000.0
# Effective bandwidth of the B-band filter in Angstroms
delta_lambda_B = 980.0
# Overall system efficiency (atmosphere * optics * detector QE)
eta = 0.4

print("This script calculates the exposure time based on the following parameters:")
print(f"Telescope Diameter: {D_m} m")
print(f"Star Apparent Magnitude (B-band): {m_B}")
print(f"Desired Magnitude Accuracy: +/- {delta_m}")
print(f"Assumed B-band Zero-Point Flux Density: {f0_lambda} photons/s/cm^2/A")
print(f"Assumed B-band Filter Bandwidth: {delta_lambda_B} A")
print(f"Assumed Overall System Efficiency: {eta}\n")


# Step 1: Calculate the total number of photons (N) required for the desired accuracy.
# The formula is N = (1.086 / delta_m)^2
print("Step 1: Calculate the required number of photons (N).")
N = (1.086 / delta_m)**2
print(f"N = (1.086 / {delta_m})^2 = {N:.2f} photons")
print("-" * 30)

# Step 2: Calculate the photon detection rate (R) from the star.
print("Step 2: Calculate the photon detection rate (R).")
# Telescope radius in cm
r_cm = (D_m * 100) / 2.0
# Telescope collecting area in cm^2
A_cm2 = math.pi * r_cm**2
print(f"Telescope Area (A) = pi * ({r_cm:.0f} cm)^2 = {A_cm2:.2f} cm^2")

# Photon rate for a magnitude 0 star (R0)
R0 = f0_lambda * A_cm2 * delta_lambda_B * eta
print(f"Rate for m_B=0 star (R0) = {f0_lambda} * {A_cm2:.2f} * {delta_lambda_B} * {eta} = {R0:.2e} photons/s")

# Photon rate for the target star (R) using the magnitude relation: R = R0 * 10^(-0.4 * m_B)
R = R0 * 10**(-0.4 * m_B)
print(f"Rate for m_B={m_B} star (R) = {R0:.2e} * 10^(-0.4 * {m_B}) = {R:.2f} photons/s")
print("-" * 30)

# Step 3: Calculate the required exposure time (t).
print("Step 3: Calculate the required exposure time (t).")
# t = N / R
t_sec = N / R
print(f"t = N / R = {N:.2f} photons / {R:.2f} photons/s = {t_sec:.2f} seconds")
print("-" * 30)

# Step 4: Round the result to the nearest integer.
final_t_sec = round(t_sec)
print(f"The final required exposure time, rounded to the nearest integer, is {final_t_sec} seconds.")
