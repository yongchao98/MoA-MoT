import math

# This script calculates the exposure time needed to achieve a specific photometric accuracy.

# --- Step 0: Define constants and problem parameters ---
mag_target = 20.0        # B-band magnitude of the star
mag_error = 0.01         # Desired magnitude accuracy in +-
telescope_diameter = 1.0 # meters, as given

# This is a standard calibrated value for the B-band zero-point flux,
# representing the number of photons per second per square meter for a magnitude 0 star,
# assuming 100% efficiency and no atmospheric extinction.
F0_B = 9.8e9

# --- Step 1: Calculate the required Signal-to-Noise Ratio (SNR) ---
# The relationship between magnitude error (dm) and SNR is SNR ≈ 1.086 / dm.
# A more precise formula is SNR = (2.5 / ln(10)) / dm.
snr_required = (2.5 / math.log(10)) / mag_error

# --- Step 2: Calculate the total number of photons (N) required ---
# For photon-counting statistics (a Poisson process), SNR = sqrt(N).
# Therefore, the total number of photons needed is N = SNR^2.
photons_needed = snr_required**2

# --- Step 3: Calculate the photon flux from the target star ---
# The flux (F) of a star with magnitude `m` is F = F0 * 10^(-m/2.5).
flux_star = F0_B * 10**(-mag_target / 2.5)

# --- Step 4: Calculate the telescope's light-collecting area ---
# The area A of a circular aperture is pi * r^2.
telescope_radius = telescope_diameter / 2.0
telescope_area = math.pi * telescope_radius**2

# --- Step 5: Calculate the required exposure time ---
# The total number of collected photons is N = flux * area * time.
# Therefore, the exposure time is t = N / (flux * area).
exposure_time = photons_needed / (flux_star * telescope_area)
rounded_time = round(exposure_time)

# --- Step 6: Display the final equation and the result ---
# Displaying the final equation with the calculated numbers
print("The exposure time (t) is calculated using the formula: t = N / (F * A)")
print("where:")
print(f" N = Total photons required to achieve the desired SNR ≈ {photons_needed:.0f}")
print(f" F = Photon flux from the m={mag_target} star ≈ {flux_star:.0f} photons/s/m²")
print(f" A = Telescope collecting area ≈ {telescope_area:.4f} m²")
print("\nSubstituting the values into the equation gives:")
print(f"t = {photons_needed:.0f} / ({flux_star:.0f} * {telescope_area:.4f})")
print(f"t ≈ {exposure_time:.2f} seconds")
print(f"\nRounded to the nearest integer, the required exposure time is {rounded_time} seconds.")