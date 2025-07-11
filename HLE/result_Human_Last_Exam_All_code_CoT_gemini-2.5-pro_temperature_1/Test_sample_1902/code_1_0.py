import math

def calculate_cumulative_dose():
    """
    Calculates the cumulative surface dose to a tissue based on ionization chamber measurements.
    """
    # --- Given Information ---
    # Convert all inputs to base SI units (meters, kilograms, seconds, etc.)

    # Ionization chamber current in Amperes (C/s)
    # 2.0 pA = 2.0 * 10^-12 A
    current_A = 2.0e-12

    # Energy to create an ion pair in air (W/e), in Joules per Coulomb
    W_air_J_per_C = 33.97

    # Density of air in kg/m^3
    # 1.293 mg/cm^3 = 1.293 * (10^-6 kg) / (10^-2 m)^3 = 1.293 kg/m^3
    density_air_kg_per_m3 = 1.293

    # Photon beam width at focus in meters
    # 0.3 mm = 0.0003 m
    beam_width_m = 0.3e-3

    # Photon beam height at focus in meters
    # 6 mm = 0.006 m
    beam_height_m = 6.0e-3

    # Length of ionization chamber in meters
    # 15.1 cm = 0.151 m
    chamber_length_m = 15.1e-2

    # Effective exposure time for a point on the sample in seconds
    exposure_time_s = 0.02

    # --- Step 1: Calculate the mass of the irradiated air ---
    beam_area_m2 = beam_width_m * beam_height_m
    irradiated_volume_m3 = beam_area_m2 * chamber_length_m
    irradiated_mass_kg = density_air_kg_per_m3 * irradiated_volume_m3

    # --- Step 2: Calculate the energy deposited per second ---
    energy_rate_J_per_s = current_A * W_air_J_per_C

    # --- Step 3: Calculate the dose rate in Gy/s (J/kg/s) ---
    dose_rate_Gy_per_s = energy_rate_J_per_s / irradiated_mass_kg

    # --- Step 4: Calculate the cumulative dose in Gy ---
    cumulative_dose_Gy = dose_rate_Gy_per_s * exposure_time_s

    # --- Print the explanation and the final equation ---
    print("The cumulative surface dose can be calculated using the formula:")
    print("Dose = ( (Current * W_air) / (Density * Beam Area * Chamber Length) ) * Exposure Time\n")
    
    print("Plugging in the values in SI units:")
    print(f"Energy Rate (J/s) = Ionization Current (A) * W_air (J/C)")
    print(f"Energy Rate (J/s) = {current_A:.1e} A * {W_air_J_per_C} J/C = {energy_rate_J_per_s:.4e} J/s\n")
    
    print(f"Irradiated Mass (kg) = Density (kg/m^3) * Beam Width (m) * Beam Height (m) * Chamber Length (m)")
    print(f"Irradiated Mass (kg) = {density_air_kg_per_m3} kg/m^3 * {beam_width_m} m * {beam_height_m} m * {chamber_length_m} m = {irradiated_mass_kg:.4e} kg\n")

    print(f"Dose Rate (Gy/s) = Energy Rate (J/s) / Irradiated Mass (kg)")
    print(f"Dose Rate (Gy/s) = {energy_rate_J_per_s:.4e} J/s / {irradiated_mass_kg:.4e} kg = {dose_rate_Gy_per_s:.4e} Gy/s\n")

    print(f"Cumulative Dose (Gy) = Dose Rate (Gy/s) * Exposure Time (s)")
    print(f"Cumulative Dose (Gy) = {dose_rate_Gy_per_s:.4e} Gy/s * {exposure_time_s} s\n")
    
    print(f"Final Cumulative Dose = {cumulative_dose_Gy:.4e} Gy")
    # Also express in more common units like milliGray (mGy) or microGray (µGy)
    print(f"Which is equal to {cumulative_dose_Gy * 1e3:.4f} mGy or {cumulative_dose_Gy * 1e6:.4f} µGy.")
    
    # Return the final value in the required format
    return cumulative_dose_Gy

# Execute the function and store the result
final_dose = calculate_cumulative_dose()
# print(f"\n<<<{final_dose:.3e}>>>") # Let's provide 3 significant figures
print(f"\n<<<{final_dose}>>>")
