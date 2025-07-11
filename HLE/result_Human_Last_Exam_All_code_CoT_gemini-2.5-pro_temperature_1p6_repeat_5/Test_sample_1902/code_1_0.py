import math

def calculate_surface_dose():
    """
    Calculates the cumulative surface dose based on synchrotron imaging parameters.
    """
    # --- Input Parameters ---
    # Beam dimensions at focus
    beam_width_focus_mm = 0.3
    beam_height_focus_mm = 6.0
    # Ionization chamber parameters
    chamber_length_cm = 15.1
    current_pA = 2.0
    # Physical properties of air
    air_density_mg_cm3 = 1.293
    # Scan parameter
    exposure_time_s = 0.02
    # Physical constant: Average energy to create an ion pair in air (in Joules/Coulomb)
    W_air_J_per_C = 33.97

    print("--- Step 1: Calculate Irradiated Air Volume ---")
    # Convert units to cm
    beam_width_focus_cm = beam_width_focus_mm / 10.0
    beam_height_focus_cm = beam_height_focus_mm / 10.0

    # Calculate beam area
    beam_area_cm2 = beam_width_focus_cm * beam_height_focus_cm
    print(f"Beam Area (cm^2) = Beam Width (cm) * Beam Height (cm)")
    print(f"Beam Area (cm^2) = {beam_width_focus_cm} cm * {beam_height_focus_cm} cm = {beam_area_cm2:.3f} cm^2")

    # Calculate irradiated volume
    irradiated_volume_cm3 = beam_area_cm2 * chamber_length_cm
    print(f"\nIrradiated Volume (cm^3) = Beam Area (cm^2) * Chamber Length (cm)")
    print(f"Irradiated Volume (cm^3) = {beam_area_cm2:.3f} cm^2 * {chamber_length_cm} cm = {irradiated_volume_cm3:.4f} cm^3")

    print("\n--- Step 2: Calculate Irradiated Air Mass ---")
    # Convert density from mg/cm^3 to g/cm^3
    air_density_g_cm3 = air_density_mg_cm3 / 1000.0

    # Calculate irradiated mass in grams
    irradiated_mass_g = irradiated_volume_cm3 * air_density_g_cm3
    print(f"Irradiated Mass (g) = Irradiated Volume (cm^3) * Air Density (g/cm^3)")
    print(f"Irradiated Mass (g) = {irradiated_volume_cm3:.4f} cm^3 * {air_density_g_cm3:.6f} g/cm^3 = {irradiated_mass_g:.6e} g")
    
    # Convert mass to kg for dose calculation
    irradiated_mass_kg = irradiated_mass_g / 1000.0
    print(f"Irradiated Mass (kg) = {irradiated_mass_g:.6e} g / 1000 = {irradiated_mass_kg:.6e} kg")

    print("\n--- Step 3: Calculate Energy Deposition Rate (Power) ---")
    # Convert current from picoamperes (pA) to amperes (A)
    current_A = current_pA * 1e-12
    
    # Calculate power in Joules/second (Watts)
    power_J_per_s = current_A * W_air_J_per_C
    print(f"Power (J/s) = Current (C/s) * W_air (J/C)")
    print(f"Power (J/s) = {current_A:.1e} C/s * {W_air_J_per_C} J/C = {power_J_per_s:.3e} J/s")

    print("\n--- Step 4: Calculate Dose Rate ---")
    # Dose Rate in Gray/second (Gy/s), where 1 Gy = 1 J/kg
    dose_rate_Gy_per_s = power_J_per_s / irradiated_mass_kg
    print(f"Dose Rate (Gy/s) = Power (J/s) / Irradiated Mass (kg)")
    print(f"Dose Rate (Gy/s) = {power_J_per_s:.3e} J/s / {irradiated_mass_kg:.6e} kg = {dose_rate_Gy_per_s:.6e} Gy/s")

    print("\n--- Step 5: Calculate Cumulative Dose ---")
    # Cumulative Dose in Gray (Gy)
    cumulative_dose_Gy = dose_rate_Gy_per_s * exposure_time_s
    print(f"Cumulative Dose (Gy) = Dose Rate (Gy/s) * Exposure Time (s)")
    print(f"Cumulative Dose (Gy) = {dose_rate_Gy_per_s:.6e} Gy/s * {exposure_time_s} s = {cumulative_dose_Gy:.6e} Gy")

    # Final result presentation
    cumulative_dose_uGy = cumulative_dose_Gy * 1e6
    print("\n--- Final Answer ---")
    print(f"The cumulative surface dose is {cumulative_dose_Gy:.3e} Gy, which is equivalent to {cumulative_dose_uGy:.2f} microGray (uGy).")

if __name__ == '__main__':
    calculate_surface_dose()
    # To conform to the output format, the final calculated value in Gy is also printed here.
    # Re-calculating to get the raw number for the final tag
    w_cm=0.3/10.0
    h_cm=6.0/10.0
    l_cm=15.1
    vol_cm3 = (w_cm*h_cm)*l_cm
    rho_g_cm3 = 1.293/1000.0
    mass_g = vol_cm3*rho_g_cm3
    mass_kg = mass_g/1000.0
    current_A = 2.0e-12
    W_air=33.97
    power_J_s = current_A * W_air
    dose_rate_Gy_s = power_J_s/mass_kg
    dose_Gy = dose_rate_Gy_s * 0.02
    print(f"\n<<<{dose_Gy}>>>")
