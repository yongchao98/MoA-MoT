import math

def calculate_exposure_time():
    """
    Calculates the exposure time needed to measure a star's B-band magnitude
    to a specified accuracy.
    """
    # ------------------
    # 1. Define Constants
    # ------------------
    delta_m = 0.01  # Desired magnitude accuracy
    m_B = 20.0      # Star B-band magnitude
    D_tel_m = 1.0   # Telescope diameter in meters

    # B-band filter characteristics and physical constants from astronomical sources
    # Using values from Bessell, M. S. (2005), ARA&A, 43, 293
    F_lambda_0 = 4.06e-9   # Zero-point flux density for Vega mag=0 (erg/s/cm^2/Angstrom)
    lambda_B_A = 4361      # Effective wavelength of B-band (Angstrom)
    delta_lambda_B_A = 970 # Effective bandwidth of B-band (Angstrom)
    h_erg_s = 6.626e-27    # Planck's constant (erg*s)
    c_cm_s = 3.0e10        # Speed of light (cm/s)

    # -----------------------------------------------------------------
    # 2. Calculate required Signal-to-Noise Ratio (S/N) and Photons (N)
    # -----------------------------------------------------------------
    # The relation between magnitude error (dm) and S/N is: dm ~ 1.0857 / (S/N)
    # For photon counting, S/N = sqrt(N), so N = (S/N)^2
    S_N_ratio = (2.5 / math.log(10)) / delta_m
    N_photons = S_N_ratio**2

    # ----------------------------------------
    # 3. Calculate Photon Collection Rate (R_p)
    # ----------------------------------------
    # a. Energy of a single B-band photon (E = hc/lambda)
    lambda_B_cm = lambda_B_A * 1e-8 # convert Angstroms to cm
    E_photon_erg = (h_erg_s * c_cm_s) / lambda_B_cm

    # b. Photon flux for a magnitude 0 star
    # Convert energy flux density to photon flux density
    f_ph_0 = F_lambda_0 / E_photon_erg  # units: photons/s/cm^2/Angstrom
    # Multiply by filter bandwidth to get total photon flux
    F_ph_0 = f_ph_0 * delta_lambda_B_A   # units: photons/s/cm^2

    # c. Photon flux for the m_B=20 star (F_m = F_0 * 10^(-m/2.5))
    F_ph_20 = F_ph_0 * 10**(-m_B / 2.5)

    # d. Telescope collecting area (A = pi * r^2)
    R_tel_cm = (D_tel_m * 100) / 2
    A_tel_cm2 = math.pi * R_tel_cm**2

    # e. Rate of photons collected (assuming 100% efficiency)
    R_p = F_ph_20 * A_tel_cm2 # units: photons/s

    # ---------------------------------
    # 4. Calculate Exposure Time (t_exp)
    # ---------------------------------
    t_exp_s = N_photons / R_p
    final_t_exp = round(t_exp_s)

    # ----------------
    # 5. Print Results
    # ----------------
    print("This script calculates the exposure time for a given star magnitude and desired precision.\n")
    print("--- Calculation Steps ---")
    print(f"1. For a magnitude error of {delta_m}, the required Signal-to-Noise ratio is {S_N_ratio:.2f}.")
    print(f"2. To achieve this S/N, the total number of photons needed (N) is {N_photons:.1f}.")
    print(f"3. The photon collection rate (R_p) for a m_B={m_B} star with a {D_tel_m}m telescope is {R_p:.2f} photons/s.")
    print("\n--- Final Equation ---")
    print(f"Exposure Time = Total Photons (N) / Photon Rate (R_p)")
    print(f"              = {N_photons:.1f} / {R_p:.2f}")
    print(f"              = {t_exp_s:.2f} seconds")
    print("\n--- Answer ---")
    print(f"The required exposure time, rounded to the nearest integer, is {final_t_exp} seconds.")


# Execute the function
calculate_exposure_time()