import numpy as np

def calculate_exposure_time():
    """
    Calculates the exposure time needed to measure a star's magnitude to a
    given precision.
    """
    # --- Step 1: Define physical and observational parameters ---
    # Desired magnitude accuracy
    delta_m = 0.01

    # Apparent magnitude of the star in the B-band
    mag_star = 20.0

    # Telescope diameter in meters
    D_tel_m = 1.0

    # Photon flux zero-point for B-band (photons/s/cm^2/Angstrom).
    # This is a standard approximate value for a magnitude 0 star.
    F0_lambda = 1000.0

    # Effective bandwidth of the B-band filter in Angstroms
    delta_lambda_B = 980.0

    # Overall system efficiency (atmosphere, telescope, detector)
    eta = 0.3

    # --- Step 2: Calculate the required Signal-to-Noise Ratio (SNR) ---
    # The relationship is delta_m approx 1.086 / SNR
    snr_factor = 2.5 / np.log(10)
    SNR_req = snr_factor / delta_m

    # --- Step 3: Calculate the total number of photons required ---
    # For photon-counting statistics, SNR = sqrt(N), so N = SNR^2
    N_photons_req = SNR_req**2

    # --- Step 4: Calculate the detected photon rate from the star ---
    # Telescope collecting area in cm^2
    A_tel_cm2 = np.pi * (D_tel_m * 100 / 2.0)**2

    # Photon rate from a magnitude 0 star, considering telescope area and efficiency
    R0_det = F0_lambda * delta_lambda_B * A_tel_cm2 * eta

    # Photon rate from the target star (m_B = 20) using the flux ratio
    R_star = R0_det * 10**(-0.4 * mag_star)

    # --- Step 5: Calculate the required exposure time ---
    # Exposure time (t) = Total photons (N) / Photon rate (R)
    t_exposure_s = N_photons_req / R_star
    t_exposure_rounded = int(round(t_exposure_s))
    
    # --- Step 6: Print the results including the final equation ---
    print(f"To achieve a magnitude precision of +/-{delta_m}, a total of {N_photons_req:.2f} photons are needed.")
    print(f"The detected photon rate from the m_B={mag_star:.0f} star is {R_star:.2f} photons/s.")
    print("-" * 50)
    print("The final equation for exposure time (t) is:")
    print(f"t = (Total Photons Needed) / (Photon Rate)")
    print(f"t = {N_photons_req:.2f} / {R_star:.2f}")
    print(f"t = {t_exposure_s:.2f} seconds")
    print("-" * 50)
    print(f"The required exposure time rounded to the nearest integer is: {t_exposure_rounded} seconds.")


if __name__ == '__main__':
    calculate_exposure_time()