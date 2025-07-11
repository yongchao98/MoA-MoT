import numpy as np

def calculate_exposure_time():
    """
    Calculates the exposure time needed to achieve a given photometric accuracy
    for a star of a certain magnitude with a given telescope.
    """

    # --- Step 0: Define Constants and Assumptions ---

    # Star and Observation Parameters
    m_B = 20.0  # Apparent B-band magnitude of the star
    delta_m = 0.01  # Desired magnitude accuracy

    # Telescope Parameters
    D = 1.0  # Telescope diameter in meters

    # Photometric and System Parameters
    # Zero-point flux for B-band: Flux from an m_B=0 star.
    # A standard value is ~1.43e6 photons/s/cm^2 (integrated over the band).
    # We convert it to photons/s/m^2 by multiplying by 10^4.
    F0_photons_per_s_per_m2 = 1.43e10

    # Assumed total system efficiency (atmosphere, optics, detector)
    eta = 0.4

    print("--- Calculation Steps ---")
    print(f"Assuming Apparent Magnitude m_B = {m_B}, Desired Accuracy delta_m = {delta_m}")
    print(f"Assuming Telescope Diameter D = {D} m, Total Efficiency eta = {eta:.2f}\n")

    # --- Step 1: Calculate the required number of photons (N) ---
    # The relation is delta_m ≈ 1.0857 / SNR. For photon counting, SNR = sqrt(N).
    # So, N = (1.0857 / delta_m)^2.
    snr_required = 1.0857 / delta_m
    N_required = snr_required**2
    print(f"Step 1: The required number of photons (N) is calculated as:")
    print(f"N = (1.0857 / {delta_m:.2f})**2 = {N_required:.1f} photons")
    print("-" * 25)

    # --- Step 2: Calculate the photon flux from the star ---
    # The flux from the star is F_star = F0 * 10^(-0.4 * m_B).
    F_star_photons_per_s_per_m2 = F0_photons_per_s_per_m2 * 10**(-0.4 * m_B)
    print("Step 2: The incident photon flux from the star (F_star) is:")
    print(f"F_star = ({F0_photons_per_s_per_m2:.2e}) * 10^(-0.4 * {m_B}) = {F_star_photons_per_s_per_m2:.2f} photons/s/m^2")
    print("-" * 25)

    # --- Step 3: Calculate the telescope collecting area (A) ---
    A = np.pi * (D / 2)**2
    print("Step 3: The telescope's collecting area (A) is:")
    print(f"A = pi * ({D} / 2)^2 = {A:.4f} m^2")
    print("-" * 25)

    # --- Step 4: Calculate the rate of detected photons (R) ---
    # R = F_star * A * eta
    R_detected_photons_per_s = F_star_photons_per_s_per_m2 * A * eta
    print("Step 4: The rate of detected photons (R) is:")
    print(f"R = F_star * A * eta = ({F_star_photons_per_s_per_m2:.2f}) * ({A:.4f}) * ({eta:.2f}) = {R_detected_photons_per_s:.2f} photons/s")
    print("-" * 25)

    # --- Step 5: Calculate the required exposure time (t) ---
    # t = N / R
    t_seconds = N_required / R_detected_photons_per_s
    t_rounded = int(round(t_seconds))
    print("Step 5: The required exposure time (t) is:")
    print(f"t = N / R = {N_required:.1f} / {R_detected_photons_per_s:.2f} = {t_seconds:.1f} seconds")
    print("-" * 25)

    print(f"\nThe final exposure time, rounded to the nearest integer, is {t_rounded} seconds.")

    return t_rounded

if __name__ == '__main__':
    final_answer = calculate_exposure_time()
    # The final answer is enclosed in <<<>>>
    print(f"\n<<< {final_answer} >>>")