import math

# This script calculates the required exposure time for a telescope to measure a star's
# magnitude to a specific precision.

# Plan:
# 1. Define the physical and observational constants.
# 2. Calculate the required Signal-to-Noise Ratio (SNR) from the desired magnitude precision.
# 3. Calculate the total number of photons (N) needed to achieve this SNR, based on Poisson statistics (SNR = sqrt(N)).
# 4. Calculate the photon flux from the star given its magnitude.
# 5. Calculate the photon collection rate of the telescope given its size.
# 6. Calculate the required exposure time by dividing the total photons needed by the collection rate.
# 7. Print the steps of the calculation and the final result rounded to the nearest integer.

# --- Step 1: Define constants ---
# Desired magnitude uncertainty
delta_m = 0.01
# Apparent B-band magnitude of the star. We assume M_B=20 refers to the apparent magnitude.
m_B = 20
# Telescope diameter in meters
D_tel = 1.0
# Zero-point photon flux for the B-band (photons/s/m^2 for a m_B=0 star)
phi_0 = 1.4e10

print("--- Calculation Steps ---")

# --- Step 2: Calculate required Signal-to-Noise Ratio (SNR) ---
# The relationship between magnitude error (delta_m) and SNR is delta_m ≈ (2.5/ln(10)) / SNR ≈ 1.086 / SNR.
# So, SNR ≈ (2.5/ln(10)) / delta_m.
snr_required = (2.5 / math.log(10)) / delta_m
print(f"1. For a magnitude uncertainty of {delta_m}, the required Signal-to-Noise Ratio (SNR) is: {snr_required:.2f}")

# --- Step 3: Calculate total number of photons needed (N) ---
# For photon-counting statistics (Poisson noise), SNR = sqrt(N), where N is the number of photons.
# Therefore, N = SNR^2.
N_photons = snr_required**2
print(f"2. To achieve this SNR, the total number of photons to be collected is: {N_photons:.2f}")

# --- Step 4: Calculate photon flux from the star ---
# The flux of the star (phi_star) is related to the zero-point flux (phi_0) by:
# phi_star = phi_0 * 10^(-m_B / 2.5)
flux_ratio = 10**(-m_B / 2.5)
phi_star = phi_0 * flux_ratio
print(f"3. The photon flux from a magnitude {m_B} star is {phi_0:.2e} * 10^(-{m_B}/2.5) = {phi_star:.2f} photons/s/m^2")

# --- Step 5: Calculate telescope's photon collection rate ---
# The collecting area of the telescope is A = pi * (D/2)^2.
A_tel = math.pi * (D_tel / 2)**2
# The photon collection rate is R = phi_star * A_tel.
photon_rate = phi_star * A_tel
print(f"4. A {D_tel}m telescope has a collecting area of {A_tel:.2f} m^2 and collects photons at a rate of: {photon_rate:.2f} photons/s")

# --- Step 6: Calculate the required exposure time ---
# The exposure time (t) is the total number of photons divided by the collection rate: t = N / R.
exposure_time = N_photons / photon_rate

print("\n--- Final Calculation ---")
print(f"The final equation is: Exposure Time = Total Photons / Photon Rate")
print(f"Exposure Time = {N_photons:.2f} photons / {photon_rate:.2f} photons/s")
print(f"Calculated exposure time: {exposure_time:.2f} seconds")

# --- Step 7: Round to the nearest integer ---
rounded_time = round(exposure_time)
print(f"\nThe required exposure time rounded to the nearest integer is: {rounded_time} seconds.")