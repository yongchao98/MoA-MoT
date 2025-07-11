import math

def calculate_surface_dose():
    """
    Calculates the cumulative surface dose based on ionization chamber measurements.
    """
    # --- Given Parameters ---
    # Beam and Chamber Geometry
    beam_width_h_mm = 0.3  # mm
    beam_width_v_mm = 6.0  # mm
    chamber_length_cm = 15.1  # cm

    # Electrical and Physical Properties
    chamber_current_pA = 2.0  # pA
    air_density_mg_cm3 = 1.293  # mg/cm^3
    W_air_over_e_J_C = 33.97  # J/C, average energy to create an ion pair in air

    # Timing
    exposure_time_s = 0.02  # s

    # --- Unit Conversions for Calculation ---
    # Convert all units to a consistent system (cm, g, s, C)
    beam_width_h_cm = beam_width_h_mm / 10.0
    beam_width_v_cm = beam_width_v_mm / 10.0
    chamber_current_A = chamber_current_pA * 1e-12  # from pA to A (C/s)
    air_density_g_cm3 = air_density_mg_cm3 / 1000.0  # from mg/cm^3 to g/cm^3

    # --- Step-by-Step Calculation ---

    # Step 1: Calculate the beam area at the focus in cm^2
    beam_area_cm2 = beam_width_h_cm * beam_width_v_cm
    print(f"Step 1: Calculate Beam Area")
    print(f"Beam Area = Horizontal Width * Vertical Width")
    print(f"            = {beam_width_h_cm:.3f} cm * {beam_width_v_cm:.1f} cm = {beam_area_cm2:.4f} cm^2\n")

    # Step 2: Calculate the volume of irradiated air in the ionization chamber in cm^3
    irradiated_volume_cm3 = beam_area_cm2 * chamber_length_cm
    print(f"Step 2: Calculate Irradiated Air Volume")
    print(f"Volume = Beam Area * Chamber Length")
    print(f"       = {beam_area_cm2:.4f} cm^2 * {chamber_length_cm:.1f} cm = {irradiated_volume_cm3:.4f} cm^3\n")

    # Step 3: Calculate the mass of the irradiated air in g and kg
    irradiated_mass_g = irradiated_volume_cm3 * air_density_g_cm3
    irradiated_mass_kg = irradiated_mass_g / 1000.0
    print(f"Step 3: Calculate Irradiated Air Mass")
    print(f"Mass = Volume * Density")
    print(f"     = {irradiated_volume_cm3:.4f} cm^3 * {air_density_g_cm3:.6f} g/cm^3 = {irradiated_mass_g:.4e} g\n")

    # Step 4: Calculate the rate of energy deposition in the air in J/s
    energy_deposition_rate_J_s = chamber_current_A * W_air_over_e_J_C
    print(f"Step 4: Calculate Energy Deposition Rate")
    print(f"Energy Rate = Current * (W_air/e)")
    print(f"            = {chamber_current_A:.1e} A * {W_air_over_e_J_C:.2f} J/C = {energy_deposition_rate_J_s:.3e} J/s\n")

    # Step 5: Calculate the dose rate in the air in Gy/s (J/kg/s)
    dose_rate_Gy_s = energy_deposition_rate_J_s / irradiated_mass_kg
    print(f"Step 5: Calculate Dose Rate in Air")
    print(f"Dose Rate = Energy Rate / Mass")
    print(f"          = {energy_deposition_rate_J_s:.3e} J/s / {irradiated_mass_kg:.4e} kg = {dose_rate_Gy_s:.4e} Gy/s\n")

    # Step 6: Calculate the cumulative dose in Gy
    # The problem assumes dose to tissue is the same as dose to air.
    cumulative_dose_Gy = dose_rate_Gy_s * exposure_time_s
    print(f"Step 6: Calculate Cumulative Surface Dose to Tissue")
    print(f"The dose to tissue is assumed to be equal to the dose in air.")
    print(f"Cumulative Dose = Dose Rate * Exposure Time")
    print(f"                = {dose_rate_Gy_s:.4e} Gy/s * {exposure_time_s:.2f} s = {cumulative_dose_Gy:.3e} Gy")
    
    return cumulative_dose_Gy

# Execute the calculation and print the final answer in the required format
final_dose = calculate_surface_dose()
print(f"\n<<<Dose: {final_dose:.3e} Gy>>>")
