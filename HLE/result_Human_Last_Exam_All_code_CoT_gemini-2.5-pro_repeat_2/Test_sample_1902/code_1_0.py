import math

def calculate_surface_dose():
    """
    Calculates the cumulative surface dose based on the provided parameters.
    """
    # --- Given Information ---
    # Convert all initial values to a consistent SI unit system (meters, kilograms, seconds, Coulombs, Joules)

    # Photon beam size at focus (horizontal by vertical)
    beam_width_focus_mm = 0.3  # in mm
    beam_height_mm = 6.0       # in mm

    # Energy of photon beam (not directly needed for this calculation method, but good practice to list)
    photon_energy_keV = 33.0   # in keV

    # Ionization chamber current
    current_pA = 2.0           # in picoamperes (pA)

    # Density of air in the ionization chamber
    density_air_mg_cm3 = 1.293 # in mg/cm^3

    # Length of ionization chamber
    chamber_length_cm = 15.1   # in cm

    # Total exposure time for a point on the tissue
    total_exposure_time_s = 0.02 # in seconds

    # Physical constant: Mean energy to create an ion pair in air per unit charge
    W_air_over_e_J_C = 33.97   # in J/C

    print("Step 1: Calculate the mass of irradiated air in the ionization chamber.")
    # Convert units to SI
    beam_width_focus_m = beam_width_focus_mm / 1000.0
    beam_height_m = beam_height_mm / 1000.0
    current_A = current_pA * 1e-12
    density_air_kg_m3 = density_air_mg_cm3 * 1000.0 # (1.293 mg/cm^3) * (1g/1000mg) * (1kg/1000g) * (100cm/1m)^3 = 1.293 kg/m^3
    chamber_length_m = chamber_length_cm / 100.0

    # Calculate beam area at focus
    beam_area_m2 = beam_width_focus_m * beam_height_m
    print(f"Photon beam area at focus = {beam_width_focus_m:.4f} m * {beam_height_m:.4f} m = {beam_area_m2:.2e} m^2")

    # Calculate the volume of irradiated air in the chamber
    irradiated_volume_m3 = beam_area_m2 * chamber_length_m
    print(f"Irradiated volume in chamber = {beam_area_m2:.2e} m^2 * {chamber_length_m:.4f} m = {irradiated_volume_m3:.2e} m^3")

    # Calculate the mass of the irradiated air
    mass_air_kg = irradiated_volume_m3 * density_air_kg_m3
    print(f"Mass of irradiated air = {irradiated_volume_m3:.2e} m^3 * {density_air_kg_m3:.4f} kg/m^3 = {mass_air_kg:.2e} kg\n")

    print("Step 2: Calculate the dose rate in air.")
    # Dose Rate in Gy/s (J/kg/s)
    # Dose Rate = (Current / Mass) * (W_air / e)
    dose_rate_air_Gy_s = (current_A / mass_air_kg) * W_air_over_e_J_C
    print(f"Dose Rate = ({current_A:.2e} C/s / {mass_air_kg:.2e} kg) * {W_air_over_e_J_C:.2f} J/C")
    print(f"Dose Rate in Air = {dose_rate_air_Gy_s:.4e} Gy/s\n")

    print("Step 3: Calculate the cumulative surface dose in tissue.")
    print("Assuming dose rate in tissue is equal to dose rate in air.")
    # The dose rate in tissue is assumed to be the same as in air
    dose_rate_tissue_Gy_s = dose_rate_air_Gy_s

    # Cumulative Dose = Dose Rate * Total Exposure Time
    cumulative_dose_Gy = dose_rate_tissue_Gy_s * total_exposure_time_s
    print(f"Cumulative Dose = {dose_rate_tissue_Gy_s:.4e} Gy/s * {total_exposure_time_s:.2f} s")
    print(f"Final Cumulative Dose = {cumulative_dose_Gy:.4e} Gy")

    # For better readability, also print in microgray (uGy)
    cumulative_dose_uGy = cumulative_dose_Gy * 1e6
    print(f"Final Cumulative Dose = {cumulative_dose_uGy:.3f} microgray\n")
    
    # Return the final answer in Gy for the final output format
    return cumulative_dose_Gy

# Run the calculation and store the final result
final_dose_in_Gy = calculate_surface_dose()

# The final answer in the required format
print(f"<<<{final_dose_in_Gy:.3e}>>>")