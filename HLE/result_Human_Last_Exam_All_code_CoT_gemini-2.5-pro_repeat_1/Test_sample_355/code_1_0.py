import math

def calculate_exposure_time():
    """
    Calculates the exposure time needed to measure a star's magnitude to a given accuracy.
    """
    # --- Problem constants ---
    magnitude_error = 0.01
    star_magnitude = 20
    telescope_diameter_m = 1.0

    # --- Standard astronomical values for the B-band ---
    # Photon flux for a magnitude 0 star in the B-band.
    photon_flux_zero_point = 1000  # photons / s / cm^2 / Angstrom
    # Effective bandwidth of the B-band filter.
    b_band_width_a = 940  # Angstroms

    # Step 1: Calculate the required signal-to-noise ratio (SNR)
    # The relationship is magnitude_error ≈ 1.086 / SNR
    snr_required = 1.086 / magnitude_error

    # Step 2: Calculate the total number of photons (N) needed
    # For photon counting statistics, SNR = sqrt(N), so N = SNR^2
    total_photons_needed = snr_required**2

    # Step 3: Calculate the photon collection rate (R) for the star
    # a) Telescope collecting area in cm^2
    telescope_radius_cm = (telescope_diameter_m * 100) / 2
    telescope_area_cm2 = math.pi * telescope_radius_cm**2

    # b) Photon rate (photons/s) collected for a magnitude 0 star
    flux_mag0_per_cm2 = photon_flux_zero_point * b_band_width_a
    rate_mag0 = flux_mag0_per_cm2 * telescope_area_cm2

    # c) Photon rate (photons/s) for the magnitude 20 star
    flux_ratio = 10**(-star_magnitude / 2.5)
    rate_mag20 = rate_mag0 * flux_ratio

    # Step 4: Calculate the required exposure time (t)
    # Total photons N = Rate R * time t  =>  t = N / R
    exposure_time_s = total_photons_needed / rate_mag20
    rounded_exposure_time = round(exposure_time_s)

    # --- Output the results and the final equation ---
    print(f"To achieve a magnitude error of {magnitude_error}, a total of {total_photons_needed:.0f} photons must be collected.")
    print(f"The telescope collects photons from the star at a rate of {rate_mag20:.2f} photons/s.")
    print("\nThe final calculation for the exposure time is:")
    print(f"Time = (Total Photons) / (Photon Rate)")
    print(f"Time = {total_photons_needed:.2f} / {rate_mag20:.2f} = {exposure_time_s:.2f} seconds")
    print(f"\nThe required exposure time, rounded to the nearest integer, is {rounded_exposure_time} seconds.")

calculate_exposure_time()
<<<160>>>