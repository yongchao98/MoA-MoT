import math

def calculate_surface_dose():
    """
    Calculates the cumulative surface dose based on provided parameters of a synchrotron imaging system.
    """
    # Given parameters from the problem description
    beam_width_focus_mm = 0.3
    beam_height_mm = 6.0
    chamber_length_cm = 15.1
    current_pA = 2.0
    rho_air_mg_cm3 = 1.293
    
    # The parameter "Ratio of beam’s width at focus to the scan step of subject per exposure time, 0.02s"
    # is interpreted as the effective total time a single point on the tissue is irradiated.
    effective_exposure_time_s = 0.02

    # Physical constant: Average energy required to produce an ion pair in dry air (W_air) divided by the elementary charge (e).
    # This value is approximately 33.97 J/C.
    W_air_over_e_JC = 33.97

    # --- Unit Conversions for Calculation ---
    beam_width_focus_cm = beam_width_focus_mm / 10.0
    beam_height_cm = beam_height_mm / 10.0
    current_A = current_pA * 1e-12  # picoamperes to amperes
    rho_air_kg_cm3 = rho_air_mg_cm3 * 1e-6  # mg/cm^3 to kg/cm^3

    print("Calculating Cumulative Surface Dose\n" + "="*40)

    # --- Step 1: Calculate the irradiated volume and mass of air ---
    print("Step 1: Calculate the mass of the irradiated air (m_air)")
    irradiated_volume_cm3 = beam_width_focus_cm * beam_height_cm * chamber_length_cm
    irradiated_mass_kg = irradiated_volume_cm3 * rho_air_kg_cm3
    print(f"Irradiated Volume = Beam Width ({beam_width_focus_cm} cm) * Beam Height ({beam_height_cm} cm) * Chamber Length ({chamber_length_cm} cm)")
    print(f"Resulting Volume = {irradiated_volume_cm3:.4f} cm^3")
    print(f"Mass = Volume ({irradiated_volume_cm3:.4f} cm^3) * Air Density ({rho_air_kg_cm3:.4e} kg/cm^3)")
    print(f"Resulting Mass = {irradiated_mass_kg:.4e} kg")
    print("-" * 40)

    # --- Step 2: Calculate the energy deposited per second (Power) ---
    print("Step 2: Calculate the energy deposited per second (Power) in the air")
    power_Js = current_A * W_air_over_e_JC
    print("Power (J/s) = Current (C/s) * (W_air / e) (J/C)")
    print(f"Power = {current_A:.1e} C/s * {W_air_over_e_JC} J/C")
    print(f"Resulting Power = {power_Js:.4e} J/s")
    print("-" * 40)

    # --- Step 3: Calculate the dose rate in air (and tissue) ---
    print("Step 3: Calculate the dose rate")
    print("Since tissue absorption is similar to air, Dose_Rate_tissue ≈ Dose_Rate_air.")
    dose_rate_Gys = power_Js / irradiated_mass_kg
    print("Dose Rate (Gy/s) = Power (J/s) / Mass (kg)")
    print(f"Dose Rate = {power_Js:.4e} J/s / {irradiated_mass_kg:.4e} kg")
    print(f"Resulting Dose Rate = {dose_rate_Gys:.4e} Gy/s (or {dose_rate_Gys * 1000:.4f} mGy/s)")
    print("-" * 40)

    # --- Step 4: Calculate the cumulative surface dose ---
    print("Step 4: Calculate the final cumulative dose")
    cumulative_dose_Gy = dose_rate_Gys * effective_exposure_time_s
    cumulative_dose_mGy = cumulative_dose_Gy * 1000
    print("Cumulative Dose (Gy) = Dose Rate (Gy/s) * Effective Exposure Time (s)")
    print(f"Cumulative Dose = {dose_rate_Gys:.4e} Gy/s * {effective_exposure_time_s} s")
    print(f"Resulting Cumulative Dose = {cumulative_dose_Gy:.4e} Gy")
    print("\n--- Final Answer ---")
    print("Final Equation: Dose (mGy) = ( (I * W/e) / (V * ρ) ) * t * 1000")
    print(f"Dose (mGy) = ( (({current_pA}e-12 A) * {W_air_over_e_JC} J/C) / (({beam_width_focus_cm} cm * {beam_height_cm} cm * {chamber_length_cm} cm) * {rho_air_mg_cm3}e-6 kg/cm^3) ) * {effective_exposure_time_s} s * 1000 mGy/Gy")
    print(f"The cumulative surface dose to the tissue is {cumulative_dose_mGy:.5f} mGy.")

    # Return the final numerical answer as requested
    return cumulative_dose_mGy

if __name__ == '__main__':
    final_dose = calculate_surface_dose()
    # The prompt requires the answer in the format <<<answer>>>
    # We print it here separately to conform to the format.
    # The main function call prints the detailed steps.
    # We round to 3 significant figures.
    print(f"\n<<<{final_dose:.3g}>>>")
