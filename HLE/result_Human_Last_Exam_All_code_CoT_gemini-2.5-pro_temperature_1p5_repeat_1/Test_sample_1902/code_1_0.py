import math

def calculate_surface_dose():
    """
    Calculates the cumulative surface dose based on ionization chamber measurements.
    """
    # --- GIVEN PARAMETERS (converted to SI units) ---

    # Photon beam size at focus (horizontal x vertical)
    beam_width_h_focus_m = 0.3 / 1000  # 0.3 mm to m
    beam_width_v_focus_m = 6.0 / 1000  # 6.0 mm to m

    # Ionization chamber parameters
    chamber_length_m = 15.1 / 100      # 15.1 cm to m
    chamber_current_A = 2.0e-12        # 2.0 pA to A

    # Physical properties of air
    # Density of air: 1.293 mg/cm^3 is equivalent to 1.293 kg/m^3
    air_density_kg_m3 = 1.293
    # Mean energy to create an ion pair in air per unit charge (W_air/e)
    W_over_e_air_J_C = 33.97           # J/C, a standard value

    # Scanning parameter
    exposure_time_s = 0.02             # 0.02 s

    # --- CALCULATION STEPS ---

    # Step 1: Calculate the mass of the irradiated air in the chamber
    beam_area_m2 = beam_width_h_focus_m * beam_width_v_focus_m
    irradiated_volume_m3 = beam_area_m2 * chamber_length_m
    mass_air_kg = irradiated_volume_m3 * air_density_kg_m3

    # Step 2: Calculate the rate of energy deposition (Power) in Watts
    power_deposited_W = chamber_current_A * W_over_e_air_J_C

    # Step 3: Calculate the absorbed dose rate in Gray per second (Gy/s)
    # Assumes dose rate in tissue is equal to dose rate in air
    dose_rate_Gy_s = power_deposited_W / mass_air_kg

    # Step 4: Calculate the final cumulative surface dose in Gray (Gy)
    cumulative_dose_Gy = dose_rate_Gy_s * exposure_time_s
    
    # --- OUTPUT RESULTS ---

    print("This script calculates the cumulative surface dose using ionization chamber data.")
    print("-" * 60)
    print("Equation for Cumulative Dose:")
    print("Dose = ( (Current * W/e) / ((Beam_H * Beam_V * Chamber_L) * Density) ) * Time")
    print("-" * 60)
    print("Substituting the given values:")
    
    # Print the full equation with all the numbers
    full_equation_str = (
        f"Dose = ( ({chamber_current_A} A * {W_over_e_air_J_C} J/C) / "
        f"(({beam_width_h_focus_m} m * {beam_width_v_focus_m} m * {chamber_length_m} m) * {air_density_kg_m3} kg/m^3) ) * "
        f"{exposure_time_s} s"
    )
    print(full_equation_str)

    # Print intermediate steps for clarity
    print("\nIntermediate Calculation Steps:")
    print(f"1. Irradiated mass of air = {mass_air_kg:.4e} kg")
    print(f"2. Energy deposition rate = {power_deposited_W:.4e} J/s")
    print(f"3. Dose rate = {dose_rate_Gy_s:.4e} Gy/s")
    
    # Print the final result
    print("-" * 60)
    print("Final Calculation:")
    print(f"Cumulative Dose = {dose_rate_Gy_s:.4e} Gy/s * {exposure_time_s} s")
    print(f"Cumulative Dose = {cumulative_dose_Gy:.4e} Gy")
    print(f"(This is equivalent to {cumulative_dose_Gy * 1e6:.3f} microGrays)")
    
    global final_answer_value
    final_answer_value = cumulative_dose_Gy

# Global variable to hold the final answer for the '<<<' tag
final_answer_value = 0
calculate_surface_dose()