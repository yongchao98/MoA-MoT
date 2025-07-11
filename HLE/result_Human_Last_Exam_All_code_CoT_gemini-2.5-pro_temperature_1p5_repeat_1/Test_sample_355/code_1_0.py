import math

def calculate_exposure_time():
    """
    Calculates the exposure time needed to measure a star's magnitude to a given accuracy.
    """
    # --- Parameters and Constants ---
    # Desired magnitude accuracy
    delta_m = 0.01
    # Apparent B-band magnitude of the star
    m_B = 20.0
    # Telescope diameter in meters
    D_meters = 1.0

    # --- Assumptions for Calculation ---
    # Standard zero-point flux for B-band (a common textbook value)
    # Units: photons / second / cm^2 / Angstrom
    phi_0 = 1000.0
    # Bandwidth of the Johnson B-band filter in Angstroms
    delta_lambda_B = 980.0
    # Assumed overall efficiency (atmosphere * telescope optics * detector QE)
    eta = 0.3

    print("Solving for the required exposure time with the following parameters:")
    print(f"  - Target Star Magnitude (M_B): {m_B}")
    print(f"  - Desired Magnitude Accuracy (δm): ±{delta_m}")
    print(f"  - Telescope Diameter: {D_meters} m")
    print(f"  - Assumed Overall Efficiency (η): {eta*100}%\n")


    # Step 1: Calculate required Signal-to-Noise Ratio (S/N)
    # The relation is δm ≈ 2.5 / (ln(10) * S/N) ≈ 1.086 / S/N
    snr = (2.5 / math.log(10)) / delta_m
    print(f"Step 1: Required Signal-to-Noise Ratio (S/N) = (2.5 / ln(10)) / {delta_m} = {snr:.2f}")

    # Step 2: Calculate required number of photons (N)
    # For photon-counting statistics, S/N = sqrt(N), so N = S/N^2
    N_photons = snr**2
    print(f"Step 2: Required number of detected photons (N) = S/N² = {snr:.2f}² = {N_photons:.2f}")

    # Step 3: Calculate the photon detection rate (R)
    # Convert telescope diameter to cm
    D_cm = D_meters * 100
    # Telescope collecting area in cm^2
    A_cm2 = math.pi * (D_cm / 2)**2
    # Photon flux from the magnitude 20 star: Flux = Flux_0 * 10^(-m / 2.5)
    phi_star = phi_0 * (10**(-m_B / 2.5))
    # Total photon detection rate: R = Flux * Area * Bandwidth * Efficiency
    photon_rate = phi_star * A_cm2 * delta_lambda_B * eta
    print(f"Step 3: Photon detection rate (R) from the star = {photon_rate:.2f} photons/second")

    # Step 4: Calculate required exposure time (t)
    exposure_time_seconds = N_photons / photon_rate
    print(f"\nStep 4: Calculate the final exposure time (t = N / R).")
    print(f"  Exposure Time = {N_photons:.2f} photons / {photon_rate:.2f} photons/s")
    print(f"  Exposure Time = {exposure_time_seconds:.2f} seconds")

    # Round the result to the nearest integer
    final_answer = round(exposure_time_seconds)
    print(f"\nThe calculated exposure time rounded to the nearest integer is {final_answer} seconds.")
    return final_answer

if __name__ == '__main__':
    result = calculate_exposure_time()
    print(f"<<<{result}>>>")