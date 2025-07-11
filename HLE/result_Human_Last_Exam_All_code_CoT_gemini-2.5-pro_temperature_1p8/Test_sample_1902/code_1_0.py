import math

def calculate_surface_dose():
    """
    Calculates the cumulative surface dose based on synchrotron imaging system parameters.
    """

    # --- Given Information and Constants ---

    # Ionization Chamber Measurements and System Parameters
    chamber_current_pA = 2.0   # pA (picoamperes)
    air_density_mg_cm3 = 1.293 # mg/cm^3
    chamber_length_cm = 15.1   # cm
    mass_atten_coeff_air_cm2_g = 0.328 # cm^2/g
    beam_width_focus_mm = 0.3  # mm
    beam_height_focus_mm = 6.0   # mm
    t_exposure_s = 0.02        # s (exposure time per point)

    # Physical Constants
    W_air_eV_per_ion_pair = 33.97  # eV required to create an ion pair in air
    e_C_per_ion_pair = 1.602e-19   # elementary charge in Coulombs per ion pair
    eV_to_J = 1.602e-19            # Conversion factor from electron-volts to Joules

    # --- Step 0: Unit Conversions ---
    
    # Convert current from pA to A (Coulombs/second)
    chamber_current_A = chamber_current_pA * 1e-12
    # Convert air density from mg/cm^3 to g/cm^3
    air_density_g_cm3 = air_density_mg_cm3 * 1e-3
    # Calculate beam area at focus and convert from mm^2 to cm^2
    beam_area_cm2 = (beam_width_focus_mm * 0.1) * (beam_height_focus_mm * 0.1)

    # --- Step 1: Calculate Energy Deposition Rate in Chamber (J/s) ---
    
    # Number of ion pairs created per second in the chamber
    ion_pairs_per_sec = chamber_current_A / e_C_per_ion_pair
    # Energy deposited per second in the chamber, in Joules/second
    energy_dep_rate_chamber_J_s = ion_pairs_per_sec * W_air_eV_per_ion_pair * eV_to_J

    # --- Step 2: Determine Total Beam Power (J/s) ---
    
    # Calculate the linear energy absorption coefficient of air (mu_en)
    mu_en_air_inv_cm = mass_atten_coeff_air_cm2_g * air_density_g_cm3
    # Calculate the fraction of the beam's energy absorbed by the air in the chamber
    fraction_absorbed = 1 - math.exp(-mu_en_air_inv_cm * chamber_length_cm)
    # Calculate the total power of the beam entering the chamber
    total_beam_power_J_s = energy_dep_rate_chamber_J_s / fraction_absorbed
    
    # --- Step 3: Calculate the Dose Rate (Gy/s) ---
    
    # Calculate beam intensity (Power/Area) in J/(s*cm^2)
    beam_intensity_J_s_cm2 = total_beam_power_J_s / beam_area_cm2
    # Calculate dose rate in J/(s*g), assuming tissue attenuation ~ air attenuation
    dose_rate_J_s_g = beam_intensity_J_s_cm2 * mass_atten_coeff_air_cm2_g
    # Convert dose rate from J/(s*g) to Gy/s (J/(s*kg))
    dose_rate_Gy_s = dose_rate_J_s_g * 1000.0  # 1000 g in 1 kg

    # --- Step 4: Calculate the Cumulative Dose (Gy) ---

    cumulative_dose_Gy = dose_rate_Gy_s * t_exposure_s

    # --- Output the Final Equation and Result ---

    print("The final calculation for the cumulative surface dose is:")
    print(f"Cumulative Dose = Dose Rate * Exposure Time")
    print(f"Cumulative Dose = {dose_rate_Gy_s:.4e} Gy/s * {t_exposure_s} s")
    
    # Presenting final answer in Gray and microgray for readability
    print(f"\nThe cumulative surface dose is {cumulative_dose_Gy:.3e} Gy.")
    print(f"This is equivalent to {cumulative_dose_Gy * 1e6:.2f} \u00B5Gy.")


# Run the calculation
calculate_surface_dose()
<<<3.88e-6>>>