import math

def calculate_cumulative_dose():
    """
    Calculates the cumulative surface dose to a tissue based on ionization chamber measurements.
    """
    # --- Given Information ---
    # Beam and Chamber properties
    beam_width_h_mm = 0.3  # mm
    beam_width_v_mm = 6.0  # mm
    chamber_length_cm = 15.1 # cm

    # Electrical and Physical Constants
    chamber_current_pA = 2.0  # pA
    air_density_mg_cm3 = 1.293 # mg/cm^3
    t_exposure_s = 0.02 # s, effective exposure time for a point on the surface

    # Constants
    e_charge = 1.602e-19 # Coulombs (C), elementary charge
    W_air_eV = 33.97     # eV/ion pair, average energy to create an ion pair in air
    
    # --- Unit Conversions to SI (meters, kilograms, seconds) ---
    # Convert pA to A (C/s)
    chamber_current_A = chamber_current_pA * 1e-12
    # Convert beam dimensions from mm to m
    beam_width_h_m = beam_width_h_mm / 1000
    beam_width_v_m = beam_width_v_mm / 1000
    # Convert chamber length from cm to m
    chamber_length_m = chamber_length_cm / 100
    # Convert air density from mg/cm^3 to kg/m^3
    # 1 mg/cm^3 = (1e-6 kg) / (1e-2 m)^3 = 1 kg/m^3
    air_density_kg_m3 = air_density_mg_cm3 * 1.0
    # Convert W_air from eV to Joules
    W_air_J = W_air_eV * e_charge
    
    # --- Step 1: Calculate Energy Deposited per Second (Power) ---
    power_deposited_J_s = (chamber_current_A / e_charge) * W_air_J

    # --- Step 2: Calculate the Mass of Irradiated Air ---
    beam_area_m2 = beam_width_h_m * beam_width_v_m
    irradiated_volume_m3 = beam_area_m2 * chamber_length_m
    mass_air_kg = irradiated_volume_m3 * air_density_kg_m3

    # --- Step 3 & 4: Calculate Dose Rate in Air (and Tissue) ---
    # Dose Rate in Gray per second (Gy/s), where 1 Gy = 1 J/kg
    # Dose rate in tissue is assumed to be the same as in air.
    dose_rate_Gy_s = power_deposited_J_s / mass_air_kg

    # --- Step 5: Calculate Cumulative Dose ---
    cumulative_dose_Gy = dose_rate_Gy_s * t_exposure_s

    # --- Output the results ---
    print("--- Calculation of Cumulative Surface Dose ---")
    
    print("\nStep 1: Energy Deposited per Second (P_dep)")
    print(f"P_dep = (Current / e) * W_air")
    print(f"P_dep = ({chamber_current_A:.1e} C/s / {e_charge:.3e} C) * {W_air_J:.3e} J")
    print(f"P_dep = {power_deposited_J_s:.4e} J/s")
    
    print("\nStep 2: Mass of Irradiated Air (m_air)")
    print(f"m_air = (Beam Area * Chamber Length) * Air Density")
    print(f"m_air = (({beam_width_h_m:.1e} m * {beam_width_v_m:.1e} m) * {chamber_length_m:.3f} m) * {air_density_kg_m3:.3f} kg/m^3")
    print(f"m_air = {mass_air_kg:.4e} kg")

    print("\nStep 3 & 4: Dose Rate (D_rate)")
    print("Dose rate in tissue is assumed equal to dose rate in air.")
    print(f"D_rate = P_dep / m_air")
    print(f"D_rate = {power_deposited_J_s:.4e} J/s / {mass_air_kg:.4e} kg")
    print(f"D_rate = {dose_rate_Gy_s:.4e} Gy/s")

    print("\nStep 5: Cumulative Dose (D_cumulative)")
    print(f"D_cumulative = D_rate * Exposure Time")
    print(f"D_cumulative = {dose_rate_Gy_s:.4e} Gy/s * {t_exposure_s} s")
    print(f"D_cumulative = {cumulative_dose_Gy:.4e} Gy")
    
    # Print the full equation with all numbers
    print("\n--- Final Equation ---")
    final_equation_str = (
        f"Dose = ( (({chamber_current_A:.1e} C/s / {e_charge:.3e} C) * {W_air_J:.3e} J) / "
        f"(({beam_width_h_m:.1e} m * {beam_width_v_m:.1e} m) * {chamber_length_m:.3f} m * {air_density_kg_m3:.3f} kg/m^3) ) * "
        f"{t_exposure_s} s"
    )
    print(final_equation_str)
    
    return cumulative_dose_Gy

# Run the calculation and store the final answer
final_answer = calculate_cumulative_dose()

# The final result in the requested format
print(f"\nFinal Answer (in Gray): {final_answer}")
print(f"<<<{final_answer:.3e}>>>")
