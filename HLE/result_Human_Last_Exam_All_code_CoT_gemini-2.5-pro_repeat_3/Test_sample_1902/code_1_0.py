import math

def calculate_cumulative_dose():
    """
    Calculates the cumulative surface dose to a tissue based on ionization chamber measurements.
    """
    # --- 1. Define and print the given parameters ---
    # Convert all initial values to a consistent unit system (cm, g, s, C, J)
    beam_width_h_mm = 0.3
    beam_width_v_mm = 6.0
    chamber_length_cm = 15.1
    current_pA = 2.0
    air_density_mg_cm3 = 1.293
    exposure_time_s = 0.02
    
    # Physical constant: Mean energy to create an ion pair in air (W/e)
    W_air_J_per_C = 33.97  # J/C

    # --- Unit Conversions ---
    beam_width_h_cm = beam_width_h_mm / 10.0
    beam_width_v_cm = beam_width_v_mm / 10.0
    current_A = current_pA * 1e-12
    air_density_g_cm3 = air_density_mg_cm3 / 1000.0

    print("--- Step 1: Calculate the mass of the irradiated air ---")
    
    # --- 2. Calculate Irradiated Volume and Mass ---
    beam_area_cm2 = beam_width_h_cm * beam_width_v_cm
    print(f"Beam cross-sectional area = {beam_width_h_cm} cm * {beam_width_v_cm} cm = {beam_area_cm2:.4f} cm^2")

    irradiated_volume_cm3 = beam_area_cm2 * chamber_length_cm
    print(f"Irradiated volume of air = {beam_area_cm2:.4f} cm^2 * {chamber_length_cm} cm = {irradiated_volume_cm3:.4f} cm^3")

    irradiated_mass_g = irradiated_volume_cm3 * air_density_g_cm3
    # Convert mass to kilograms for dose calculation in Gray (J/kg)
    irradiated_mass_kg = irradiated_mass_g / 1000.0
    print(f"Mass of irradiated air = {irradiated_volume_cm3:.4f} cm^3 * {air_density_g_cm3:.6f} g/cm^3 = {irradiated_mass_kg:.4e} kg")

    print("\n--- Step 2: Calculate the dose rate ---")
    
    # --- 3. Calculate Energy Deposition Rate (Power) ---
    power_deposited_J_s = current_A * W_air_J_per_C
    print(f"Energy deposition rate (Power) = Ionization Current * (W/e)_air")
    print(f"Power = {current_A:.2e} A * {W_air_J_per_C} J/C = {power_deposited_J_s:.4e} J/s")
    
    # --- 4. Calculate Dose Rate ---
    # Dose Rate = Power / Mass. This is the dose rate in air.
    # We assume Dose_Rate_tissue is approximately Dose_Rate_air.
    dose_rate_Gy_s = power_deposited_J_s / irradiated_mass_kg
    print(f"Dose Rate = Power / Mass = {power_deposited_J_s:.4e} J/s / {irradiated_mass_kg:.4e} kg = {dose_rate_Gy_s:.4e} Gy/s")
    
    print("\n--- Step 3: Calculate the cumulative dose ---")
    
    # --- 5. Calculate Cumulative Dose ---
    cumulative_dose_Gy = dose_rate_Gy_s * exposure_time_s
    print(f"The cumulative dose is the dose rate multiplied by the exposure time.")
    print(f"Final Cumulative Dose = {dose_rate_Gy_s:.4e} Gy/s * {exposure_time_s} s = {cumulative_dose_Gy:.4e} Gy")
    
    # Return the final numerical answer for the submission format
    return cumulative_dose_Gy

# Execute the calculation and print the final answer in the required format
final_dose = calculate_cumulative_dose()
# The final answer is typically requested in a specific format.
# The calculation results in ~3.87e-6 Gy, which is 3.87 microGray.
# Printing the numerical value as the final line.
print(f"\n<<<The cumulative surface dose is {final_dose:.3e} Gy>>>")

# The problem asks for the answer in format <<<answer>>>. The final answer is the numerical value.
# For example, if the value is 3.865e-6, the answer is 3.865e-6.
final_answer_value = f"{final_dose:.4e}"
# Let's reformat to a more standard scientific notation if needed, but this is fine.
# Let's be more precise
final_answer_value = "{:.4e}".format(final_dose)
# The prompt is a bit ambiguous, let's provide the number directly.
# <<<3.8653e-06>>>
