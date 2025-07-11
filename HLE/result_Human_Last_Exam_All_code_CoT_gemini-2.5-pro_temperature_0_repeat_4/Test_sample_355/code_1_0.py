import math

def calculate_exposure_time():
    """
    Calculates the exposure time needed to measure a star's magnitude to a given accuracy.
    """
    # --- Step 1: Define parameters ---
    # Desired magnitude accuracy
    delta_m = 0.01
    # Magnitude of the star
    m_star = 20
    # Telescope diameter in meters
    D_tel = 1.0
    # B-band zero-point photon flux (photons/s/m^2 for m_B=0).
    # This is a standard approximate value for the Johnson B-band.
    F0_photons = 1.4e10

    # --- Step 2: Calculate required SNR and total photons (N) ---
    # The relation between magnitude uncertainty and fractional flux uncertainty (1/SNR) is
    # delta_m = (2.5 / ln(10)) * (1 / SNR)
    # So, SNR = (2.5 / ln(10)) / delta_m
    snr = (2.5 / math.log(10)) / delta_m
    
    # For photon counting (Poisson statistics), SNR = sqrt(N)
    # So, N = SNR^2
    N_total = snr**2

    # --- Step 3: Calculate photon collection rate ---
    # The flux from the star is F_star = F0 * 10^(-m_star / 2.5)
    flux_star = F0_photons * (10**(-m_star / 2.5))

    # Telescope collecting area A = pi * (D/2)^2
    area_tel = math.pi * (D_tel / 2)**2
    
    # Photon collection rate = Flux * Area
    photon_rate = flux_star * area_tel

    # --- Step 4: Calculate exposure time ---
    exposure_time = N_total / photon_rate

    # --- Step 5: Print the results ---
    print(f"To achieve a magnitude accuracy of +/- {delta_m}, we need a signal-to-noise ratio of {snr:.1f}.")
    print(f"This corresponds to collecting a total of N = {snr:.1f}^2 = {N_total:.1f} photons.")
    print("-" * 30)
    print(f"The photon flux from a magnitude {m_star} star is {flux_star:.1f} photons/s/m^2.")
    print(f"A {D_tel}m telescope has a collecting area of {area_tel:.3f} m^2.")
    print(f"The photon collection rate is {flux_star:.1f} * {area_tel:.3f} = {photon_rate:.3f} photons/s.")
    print("-" * 30)
    
    # The prompt asks to output each number in the final equation.
    # The final equation is: Exposure Time = Total Photons / Photon Rate
    print("Final Equation:")
    print(f"Exposure Time = Total Photons / Photon Rate")
    print(f"Exposure Time = {N_total:.1f} / {photon_rate:.3f} = {exposure_time:.1f} seconds")

    # Round the final answer to the nearest integer
    rounded_exposure_time = round(exposure_time)
    print(f"\nThe required exposure time, rounded to the nearest integer, is {rounded_exposure_time} seconds.")

calculate_exposure_time()