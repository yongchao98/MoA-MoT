import math

def calculate_surface_dose():
    """
    Calculates the cumulative surface dose based on ionization chamber measurements.
    """
    # --- 1. Define Constants and Inputs ---
    # Photon beam size at focus (horizontal x vertical)
    beam_width_h_mm = 0.3  # mm
    beam_width_v_mm = 6.0  # mm
    
    # Ionization chamber properties
    chamber_length_cm = 15.1  # cm
    current_pA = 2.0  # pA
    
    # Air properties
    density_air_mg_cm3 = 1.293  # mg/cm^3
    
    # Standard physical constant: average energy to create an ion pair in dry air
    W_air_per_e_J_C = 33.97  # J/C
    
    # The cumulative exposure time for a point on the surface, interpreted from the problem statement.
    exposure_time_s = 0.02  # s

    # --- 2. Unit Conversions for Consistency ---
    # Convert all length units to cm, mass to g, current to A
    beam_width_h_cm = beam_width_h_mm / 10.0
    beam_width_v_cm = beam_width_v_mm / 10.0
    current_A = current_pA * 1e-12  # picoamperes to amperes
    density_air_g_cm3 = density_air_mg_cm3 / 1000.0 # mg/cm^3 to g/cm^3

    # --- 3. Calculate Mass of Irradiated Air ---
    # The mass of the air in the chamber's active volume that is hit by the beam.
    beam_area_cm2 = beam_width_h_cm * beam_width_v_cm
    irradiated_volume_cm3 = beam_area_cm2 * chamber_length_cm
    mass_air_g = irradiated_volume_cm3 * density_air_g_cm3
    mass_air_kg = mass_air_g / 1000.0  # Convert mass to kilograms for dose calculation

    # --- 4. Calculate Total Energy Deposited in Air ---
    # Total charge collected by the chamber during the exposure.
    total_charge_C = current_A * exposure_time_s
    
    # Total energy deposited in the air, creating the measured charge.
    energy_deposited_J = total_charge_C * W_air_per_e_J_C

    # --- 5. Calculate the Cumulative Surface Dose ---
    # Dose = Energy / Mass. Unit: J/kg = Gray (Gy).
    # Assuming Dose_tissue is equivalent to Dose_air.
    if mass_air_kg == 0:
        dose_Gy = 0
    else:
        dose_Gy = energy_deposited_J / mass_air_kg
    
    # Convert dose to a more common unit, milliGray (mGy).
    dose_mGy = dose_Gy * 1000

    # --- 6. Print the Detailed Calculation and Result ---
    print("Step 1: Calculate the mass of irradiated air (m_air)")
    print(f"  - Beam area = {beam_width_h_cm:.2f} cm * {beam_width_v_cm:.1f} cm = {beam_area_cm2:.4f} cm^2")
    print(f"  - Irradiated volume = {beam_area_cm2:.4f} cm^2 * {chamber_length_cm} cm = {irradiated_volume_cm3:.4f} cm^3")
    print(f"  - Mass of air = {irradiated_volume_cm3:.4f} cm^3 * {density_air_g_cm3:.6f} g/cm^3 = {mass_air_g:.4e} g")
    print(f"  - m_air = {mass_air_kg:.4e} kg\n")

    print("Step 2: Calculate the energy deposited (E_dep)")
    print(f"  - Total charge = {current_A:.1e} C/s * {exposure_time_s} s = {total_charge_C:.1e} C")
    print(f"  - Energy deposited = {total_charge_C:.1e} C * {W_air_per_e_J_C} J/C = {energy_deposited_J:.4e} J\n")

    print("Step 3: Calculate the cumulative surface dose (Dose = E_dep / m_air)")
    print(f"  Dose = {energy_deposited_J:.4e} J / {mass_air_kg:.4e} kg")
    print(f"  Dose = {dose_Gy:.4e} Gy")
    print(f"  Dose = {dose_mGy:.4f} mGy\n")
    
    # Return final value for the answer tag
    return dose_mGy

if __name__ == '__main__':
    final_dose_mGy = calculate_surface_dose()
    print(f"The final calculated cumulative surface dose is {final_dose_mGy:.4f} mGy.")
    print(f"<<<{final_dose_mGy:.4f}>>>")