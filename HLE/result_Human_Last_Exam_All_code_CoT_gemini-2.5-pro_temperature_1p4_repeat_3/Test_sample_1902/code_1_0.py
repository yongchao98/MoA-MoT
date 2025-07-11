import math

def calculate_cumulative_dose():
    """
    Calculates the cumulative surface dose based on synchrotron experimental parameters.
    """
    # --- Given Information ---
    # Photon beam size at focus (horizontal by vertical)
    beam_width_focus_mm = 0.3  # mm
    beam_height_mm = 6.0      # mm

    # Length of ionization chamber
    chamber_length_cm = 15.1  # cm

    # Ionization chamber current
    current_pA = 2.0  # pA

    # Density of air in the ionization chamber
    air_density_mg_cm3 = 1.293  # mg/cm^3

    # Exposure time per scan step (dwell time)
    t_dwell_s = 0.02  # s

    # --- Physical Constants ---
    # Average energy required to produce an ion pair in air (W_air/e)
    W_air_div_e_JC = 33.97  # J/C

    # --- Unit Conversions ---
    # Convert all units to a consistent system for calculation (cm, g, s)
    # and final conversion to SI for dose (Gy = J/kg).

    # Beam dimensions in cm
    beam_width_focus_cm = beam_width_focus_mm / 10.0
    beam_height_cm = beam_height_mm / 10.0

    # Current in C/s (Amperes)
    current_A = current_pA * 1e-12

    # Air density in g/cm^3
    air_density_g_cm3 = air_density_mg_cm3 / 1000.0

    # --- Step-by-Step Calculation ---

    print("--- Calculation of Cumulative Surface Dose ---")

    # 1. Calculate the mass of the irradiated air.
    beam_area_cm2 = beam_width_focus_cm * beam_height_cm
    air_volume_cm3 = beam_area_cm2 * chamber_length_cm
    air_mass_g = air_volume_cm3 * air_density_g_cm3
    air_mass_kg = air_mass_g / 1000.0
    
    print("\nStep 1: Calculate the mass of irradiated air.")
    print(f"The mass is calculated as: Volume * Density")
    print(f"Volume = (Beam Width: {beam_width_focus_cm} cm * Beam Height: {beam_height_cm} cm) * Chamber Length: {chamber_length_cm} cm = {air_volume_cm3:.4f} cm^3")
    print(f"Mass = Volume: {air_volume_cm3:.4f} cm^3 * Air Density: {air_density_g_cm3} g/cm^3 = {air_mass_kg:.4e} kg")

    # 2. Calculate the dose rate in air (which equals dose rate in tissue).
    energy_dep_rate_Js = current_A * W_air_div_e_JC
    dose_rate_Gys = energy_dep_rate_Js / air_mass_kg
    
    print("\nStep 2: Calculate the dose rate.")
    print(f"The dose rate is calculated as: (Energy Deposited per Second) / Mass")
    print(f"Energy Deposited per Second = Current: {current_A:.1e} C/s * (W_air/e): {W_air_div_e_JC} J/C = {energy_dep_rate_Js:.4e} J/s")
    print(f"Dose Rate = Energy Deposited per Second: {energy_dep_rate_Js:.4e} J/s / Mass: {air_mass_kg:.4e} kg = {dose_rate_Gys:.4f} Gy/s")

    # 3. Calculate the final cumulative dose.
    cumulative_dose_Gy = dose_rate_Gys * t_dwell_s
    
    print("\nStep 3: Calculate the cumulative dose.")
    print(f"The cumulative dose is calculated as: Dose Rate * Exposure Time")
    print(f"Cumulative Dose = {dose_rate_Gys:.4f} Gy/s * {t_dwell_s} s")

    print(f"\nFinal Result:")
    print(f"The cumulative surface dose is {cumulative_dose_Gy:.6f} Gy.")

    return cumulative_dose_Gy

# Execute the calculation and store the final answer.
final_answer = calculate_cumulative_dose()
# The final answer is requested in a specific format.
# We will output the dose in Gray (Gy).
print(f"<<<{final_answer:.4f}>>>")
