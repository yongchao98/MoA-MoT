import math

# --- Given Parameters ---
# Apparent magnitude of the star. We assume the given M_B=20 is the apparent magnitude m_B.
m_B = 20.0
# Diameter of the telescope in meters
D_telescope = 1.0
# Desired magnitude uncertainty
delta_m = 0.01

# --- Physical Constants ---
# B-band flux for a magnitude 0 star (photons / s / m^2)
# This is a standard approximate value.
F0_B = 9.4e9

# --- Step 1: Calculate required Signal-to-Noise Ratio (SNR) ---
# The relationship between magnitude uncertainty (delta_m) and fractional flux uncertainty (delta_F / F)
# is delta_m = (2.5 / ln(10)) * (delta_F / F).
# The SNR is the inverse of the fractional flux uncertainty, SNR = F / delta_F.
# So, SNR = (2.5 / ln(10)) / delta_m
snr_required = (2.5 / math.log(10)) / delta_m
print(f"To achieve a magnitude uncertainty of {delta_m}, we need a Signal-to-Noise Ratio (SNR) of: {snr_required:.2f}")

# --- Step 2: Calculate total number of photons (N) needed ---
# For photon counting (Poisson statistics), SNR = sqrt(N).
# Therefore, N = SNR^2.
N_photons = snr_required**2
print(f"The total number of photons that need to be collected is: {N_photons:.2f}")

# --- Step 3: Calculate the photon flux from the star ---
# The flux of a star of magnitude m_B is F_star = F0_B * 10^(-m_B / 2.5)
flux_star = F0_B * (10**(-m_B / 2.5))
print(f"The photon flux from a B={m_B} star is: {flux_star:.2f} photons/s/m^2")

# --- Step 4: Calculate the photon collection rate of the telescope ---
# Telescope collecting area A = pi * (D/2)^2
radius = D_telescope / 2.0
area_telescope = math.pi * (radius**2)
print(f"The collecting area of a {D_telescope}m telescope is: {area_telescope:.4f} m^2")

# Photon collection rate R = F_star * A
photon_rate = flux_star * area_telescope
print(f"The telescope collects photons from the star at a rate of: {photon_rate:.2f} photons/s")

# --- Step 5: Calculate the required exposure time ---
# Exposure time T = N_photons / photon_rate
exposure_time = N_photons / photon_rate
rounded_time = round(exposure_time)

print("\n--- Final Calculation ---")
print(f"The required exposure time is the total number of photons divided by the collection rate.")
print(f"Exposure Time = {N_photons:.0f} photons / {photon_rate:.2f} photons/s")
print(f"Exposure Time = {exposure_time:.2f} seconds")
print(f"\nRounded to the nearest integer, the required exposure time is {rounded_time} seconds.")
