import math

def calculate_cumulative_dose():
    """
    Calculates the cumulative surface dose based on provided parameters.
    """
    # --- Given parameters ---
    beam_width_h_mm = 0.3      # horizontal beam width in mm
    beam_width_v_mm = 6        # vertical beam width in mm
    ionization_current_pA = 2.0  # pA
    density_air_mg_cm3 = 1.293 # mg/cm^3
    chamber_length_cm = 15.1   # cm
    exposure_time_s = 0.02     # s

    # --- Physical constants ---
    # W_air/e: Average energy required to create an ion pair in air, divided by the elementary charge
    W_air_over_e_J_C = 33.97   # J/C

    # --- Unit Conversions ---
    # Convert beam dimensions from mm to cm
    beam_width_h_cm = beam_width_h_mm / 10.0
    beam_width_v_cm = beam_width_v_mm / 10.0

    # Convert ionization current from picoamperes (pA) to amperes (A, or C/s)
    ionization_current_A = ionization_current_pA * 1e-12

    # Convert air density from mg/cm^3 to kg/cm^3 for SI unit compatibility
    density_air_kg_cm3 = density_air_mg_cm3 * 1e-6

    # --- Step 1: Calculate the mass of the irradiated air (m) ---
    # Calculate the cross-sectional area of the beam in cm^2
    beam_area_cm2 = beam_width_h_cm * beam_width_v_cm

    # Calculate the volume of air irradiated by the beam in cm^3
    irradiated_volume_cm3 = beam_area_cm2 * chamber_length_cm

    # Calculate the mass of the irradiated air in kg
    mass_air_kg = irradiated_volume_cm3 * density_air_kg_cm3

    # --- Step 2: Calculate the dose rate (D_dot) in Gy/s ---
    # Dose rate is energy absorbed per unit mass per unit time.
    # Energy absorbed per time = Current * (W_air/e)
    # Dose rate = (Current * (W_air/e)) / mass
    dose_rate_Gy_s = (ionization_current_A * W_air_over_e_J_C) / mass_air_kg

    # --- Step 3: Calculate the cumulative dose (D) in Gy ---
    # Cumulative dose is the dose rate multiplied by the exposure time for a single point.
    cumulative_dose_Gy = dose_rate_Gy_s * exposure_time_s

    # --- Output Results ---
    print("Equation for Cumulative Dose (D):")
    print("D = [(I * (W/e)) / (A_beam * L * rho)] * t_exp")
    print("\nWhere:")
    print(f"I (Current) = {ionization_current_A:.1e} C/s")
    print(f"W/e = {W_air_over_e_J_C} J/C")
    print(f"A_beam (Beam Area) = {beam_width_h_cm} cm * {beam_width_v_cm} cm = {beam_area_cm2:.4f} cm^2")
    print(f"L (Chamber Length) = {chamber_length_cm} cm")
    print(f"rho (Air Density) = {density_air_kg_cm3:.4e} kg/cm^3")
    print(f"t_exp (Exposure Time) = {exposure_time_s} s")
    
    print("\nCalculation of Dose Rate:")
    print(f"Dose Rate = ({ionization_current_A:.1e} C/s * {W_air_over_e_J_C} J/C) / ({beam_area_cm2:.4f} cm^2 * {chamber_length_cm} cm * {density_air_kg_cm3:.4e} kg/cm^3)")
    print(f"Dose Rate = {dose_rate_Gy_s:.4e} Gy/s")
    
    print("\nCalculation of Cumulative Dose:")
    print(f"Cumulative Dose = {dose_rate_Gy_s:.4e} Gy/s * {exposure_time_s} s")

    print(f"\nThe final cumulative surface dose is: {cumulative_dose_Gy:.3e} Gy")
    # Using scientific notation for the final answer
    # To provide a single number for the final answer format
    final_answer = "{:.3e}".format(cumulative_dose_Gy)
    return final_answer

# Run the calculation and print the result
final_answer_value = calculate_cumulative_dose()
print(f"<<<{final_answer_value}>>>")
