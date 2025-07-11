import math

def calculate_surface_dose():
    """
    Calculates the cumulative surface dose based on synchrotron experimental parameters.
    """
    # --- Input Parameters and Constants ---

    # Ionization chamber current in Amperes (C/s)
    # Given: 2.0 pA = 2.0e-12 A
    ionization_current_A = 2.0e-12

    # Photon beam size at focus in cm
    # Given: 0.3 mm by 6 mm = 0.03 cm by 0.6 cm
    beam_width_h_cm = 0.03
    beam_width_v_cm = 0.6

    # Length of the ionization chamber in cm
    # Given: 15.1 cm
    chamber_length_cm = 15.1

    # Density of air in g/cm^3
    # Given: 1.293 mg/cm^3 = 1.293e-3 g/cm^3
    air_density_g_per_cm3 = 1.293e-3

    # Energy absorption mass attenuation coefficient of air in cm^2/g
    # Given: 0.328 cm^2/g
    mu_en_rho_air_cm2_per_g = 0.328

    # The effective time a single point is exposed to the beam in seconds
    # Given: Ratio interpretation gives an effective exposure time of 0.02 s
    exposure_time_s = 0.02

    # Physical constant: Mean energy to create an ion pair in air (W/e) in J/C
    # W_air = 33.97 eV/ip, e = 1.602e-19 C/ip -> W/e = 33.97 J/C
    W_air_J_per_C = 33.97

    # --- Step 1: Calculate Absorbed Power in Chamber ---
    # Power_absorbed (J/s) = Current (C/s) * Energy_per_charge (J/C)
    power_absorbed_J_s = ionization_current_A * W_air_J_per_C

    # --- Step 2: Calculate Incident Beam Power ---
    # Attenuation factor x = (μ_en/ρ) * ρ * L
    x = mu_en_rho_air_cm2_per_g * air_density_g_per_cm3 * chamber_length_cm
    # Fraction of energy absorbed, f = 1 - exp(-x)
    fraction_absorbed = 1 - math.exp(-x)
    # Incident Power = Absorbed Power / Fraction Absorbed
    power_incident_J_s = power_absorbed_J_s / fraction_absorbed

    # --- Step 3: Calculate Beam Intensity ---
    # Area = horizontal_width * vertical_width
    area_focus_cm2 = beam_width_h_cm * beam_width_v_cm
    # Intensity = Power / Area
    intensity_J_per_s_cm2 = power_incident_J_s / area_focus_cm2

    # --- Step 4: Calculate Surface Dose Rate ---
    # Dose Rate (J/s/g) = Intensity (J/s/cm^2) * Mass Attenuation Coeff (cm^2/g)
    dose_rate_J_per_s_g = intensity_J_per_s_cm2 * mu_en_rho_air_cm2_per_g
    # Convert Dose Rate to Gy/s (1 Gy = 1 J/kg, so multiply by 1000 g/kg)
    dose_rate_Gy_per_s = dose_rate_J_per_s_g * 1000

    # --- Step 5: Calculate Cumulative Dose ---
    # Cumulative Dose (Gy) = Dose Rate (Gy/s) * Exposure Time (s)
    cumulative_dose_Gy = dose_rate_Gy_per_s * exposure_time_s

    # --- Print the Final Calculation ---
    print("The final calculation is:")
    print(f"Cumulative Dose = Dose Rate x Exposure Time")
    print(f"Cumulative Dose = {dose_rate_Gy_per_s:.4e} Gy/s * {exposure_time_s} s")
    print(f"Cumulative Dose = {cumulative_dose_Gy:.4e} Gy")
    # Also printing in milliGrays for context
    print(f"This is equivalent to {cumulative_dose_Gy * 1000:.4f} mGy")
    
    # Store the final numerical answer for the designated output format
    global final_answer
    final_answer = cumulative_dose_Gy

# Run the calculation
final_answer = 0
calculate_surface_dose()