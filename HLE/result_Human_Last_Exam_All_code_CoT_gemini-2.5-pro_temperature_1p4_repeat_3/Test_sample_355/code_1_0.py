import math

def calculate_exposure_time():
    """
    Calculates the exposure time needed to measure a star's magnitude to a given accuracy.
    """
    # --- Input Parameters ---
    # The apparent magnitude of the star. We assume M_B=20 is a typo for m_B=20.
    m_B = 20.0
    # The desired magnitude accuracy.
    delta_m = 0.01
    # The diameter of the telescope in meters.
    D = 1.0
    # Zero-point flux for the B-band (a standard rule-of-thumb value).
    # Units are photons / s / m^2 for a magnitude 0 star.
    F0 = 1e10

    # --- Step 1: Calculate the required number of photons (N) ---
    # The relationship is delta_m ≈ 1.0857 / SNR, and for photon counting, SNR = sqrt(N).
    # Rearranging for N gives: N = (1.0857 / delta_m)^2
    N = (1.0857 / delta_m)**2

    print("Step 1: Calculate total photons needed (N)")
    print(f"   N = (1.0857 / {delta_m})**2")
    print(f"   N = {N:.1f} photons\n")

    # --- Step 2: Calculate the photon collection rate (R) ---
    # a) Calculate the photon flux from the star (F_star) in photons/s/m^2
    F_star = F0 * 10**(-0.4 * m_B)

    # b) Calculate the telescope's collecting area (A) in m^2
    radius = D / 2.0
    A = math.pi * radius**2

    # c) Calculate the final collection rate (R) in photons/s
    R = F_star * A

    print("Step 2: Calculate photon collection rate (R)")
    print("   Final Equation for Rate: R = F0 * 10**(-0.4 * m_B) * pi * (D/2)**2")
    print(f"   R = {F0:.0e} * 10**(-0.4 * {m_B}) * {math.pi:.4f} * ({radius})**2")
    print(f"   R = {R:.2f} photons/s\n")


    # --- Step 3: Calculate the required exposure time (t) ---
    # The exposure time is the total number of photons needed divided by the collection rate.
    t = N / R
    t_rounded = int(round(t))

    print("Step 3: Calculate required exposure time (t)")
    print("   Final Equation for Time: t = N / R")
    print(f"   t = {N:.1f} / {R:.2f}")
    print(f"   t = {t:.2f} seconds\n")

    print("--- Final Answer ---")
    print(f"The required exposure time, rounded to the nearest integer, is {t_rounded} seconds.")

calculate_exposure_time()
<<<150>>>