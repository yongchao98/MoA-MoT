import math

def calculate_surface_dose():
    """
    Calculates the cumulative surface dose based on ionization chamber measurements.
    """
    # 1. Given information and constants
    #
    # Given data from the problem description
    beam_width_focus_mm = 0.3
    beam_height_focus_mm = 6.0
    ion_current_pA = 2.0
    air_density_mg_cm3 = 1.293
    chamber_length_cm = 15.1
    exposure_time_s = 0.02
    
    # Physical constant: Average energy to create an ion pair in air (W/e) in Joules/Coulomb
    W_air_over_e_J_C = 33.97

    # 2. Perform unit conversions for calculation
    #
    # Convert current from picoamperes (pA) to amperes (A, or C/s)
    ion_current_A = ion_current_pA * 1e-12
    # Convert beam dimensions from mm to cm
    beam_width_focus_cm = beam_width_focus_mm / 10.0
    beam_height_focus_cm = beam_height_focus_mm / 10.0
    # Convert air density from mg/cm^3 to kg/cm^3
    air_density_kg_cm3 = air_density_mg_cm3 / 1000.0 / 1000.0

    # 3. Step-by-step calculation
    #
    # a) Calculate the energy deposited per second (Power) in Joules/second
    power_J_s = ion_current_A * W_air_over_e_J_C

    # b) Calculate the volume of the irradiated air in the chamber in cm^3
    beam_area_cm2 = beam_width_focus_cm * beam_height_focus_cm
    irradiated_volume_cm3 = beam_area_cm2 * chamber_length_cm

    # c) Calculate the mass of the irradiated air in kg
    mass_air_kg = air_density_kg_cm3 * irradiated_volume_cm3

    # d) Calculate the dose rate in Gray per second (Gy/s)
    dose_rate_Gy_s = power_J_s / mass_air_kg

    # e) Calculate the final cumulative dose in Gray (Gy)
    cumulative_dose_Gy = dose_rate_Gy_s * exposure_time_s

    # 4. Print the results
    #
    print("--- Input Values ---")
    print(f"Ionization Chamber Current: {ion_current_pA} pA")
    print(f"Photon Beam Size at Focus: {beam_width_focus_mm} mm x {beam_height_focus_mm} mm")
    print(f"Length of Ionization Chamber: {chamber_length_cm} cm")
    print(f"Density of Air: {air_density_mg_cm3} mg/cm^3")
    print(f"Exposure Time: {exposure_time_s} s\n")
    
    print("--- Calculation Steps ---")
    print(f"1. Dose Rate Calculation:")
    print(f"   - Energy Rate = {ion_current_A:.2e} C/s * {W_air_over_e_J_C} J/C = {power_J_s:.4e} J/s")
    print(f"   - Irradiated Mass = ({beam_width_focus_cm} cm * {beam_height_focus_cm} cm * {chamber_length_cm} cm) * {air_density_kg_cm3:.4e} kg/cm^3 = {mass_air_kg:.4e} kg")
    print(f"   - Dose Rate = {power_J_s:.4e} J/s / {mass_air_kg:.4e} kg = {dose_rate_Gy_s:.4e} Gy/s\n")
    
    print("2. Final Cumulative Dose Calculation:")
    print(f"   Cumulative Dose = Dose Rate * Exposure Time")
    print(f"   Cumulative Dose = {dose_rate_Gy_s:.4e} Gy/s * {exposure_time_s} s\n")
    
    print("--- Final Answer ---")
    print(f"The cumulative surface dose is: {cumulative_dose_Gy:.4e} Gy")


calculate_surface_dose()