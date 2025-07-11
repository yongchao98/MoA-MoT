import numpy as np

def calculate_cumulative_dose():
    """
    Calculates the cumulative surface dose to a tissue imaged using a multiple energy imaging system.
    """
    # Step 1: Define constants and given parameters in consistent SI units where applicable.
    # Beam dimensions at focus
    beam_width_focus_mm = 0.3  # mm
    beam_height_focus_mm = 6.0   # mm
    # Convert beam dimensions to meters for SI consistency later
    beam_width_focus_m = beam_width_focus_mm / 1000.0
    beam_height_focus_m = beam_height_focus_mm / 1000.0

    # Ionization chamber parameters
    chamber_current_pA = 2.0  # pA
    chamber_length_cm = 15.1  # cm
    air_density_mg_cm3 = 1.293 # mg/cm^3

    # Physical constants
    # Mean energy to create an ion pair in air (W/e) in J/C
    W_e_air_J_C = 33.97
    
    # Exposure time
    exposure_time_s = 0.02  # seconds

    # Step 2: Calculate the mass of the irradiated air
    # Convert units to a consistent system (e.g., cm and g)
    beam_width_focus_cm = beam_width_focus_mm / 10.0
    beam_height_focus_cm = beam_height_focus_mm / 10.0
    beam_area_focus_cm2 = beam_width_focus_cm * beam_height_focus_cm
    
    # Volume of air irradiated in the chamber
    irradiated_volume_cm3 = beam_area_focus_cm2 * chamber_length_cm
    
    # Mass of the irradiated air
    air_density_g_cm3 = air_density_mg_cm3 / 1000.0 # from mg/cm^3 to g/cm^3
    irradiated_mass_g = irradiated_volume_cm3 * air_density_g_cm3
    irradiated_mass_kg = irradiated_mass_g / 1000.0 # from g to kg

    # Step 3: Calculate the dose rate in air (and tissue)
    # Convert current from pA to A (C/s)
    chamber_current_A = chamber_current_pA * 1e-12

    # Rate of energy absorption (J/s)
    energy_absorption_rate_J_s = chamber_current_A * W_e_air_J_C
    
    # Dose rate (Gy/s = J/kg/s)
    dose_rate_Gy_s = energy_absorption_rate_J_s / irradiated_mass_kg

    # Step 4: Calculate the cumulative dose
    # The dose rate in tissue is assumed to be the same as in air.
    # Cumulative Dose = Dose Rate * Exposure Time
    cumulative_dose_Gy = dose_rate_Gy_s * exposure_time_s

    # Print the final calculation as an equation
    print("The cumulative surface dose is calculated by multiplying the dose rate by the exposure time.")
    print("Final Equation:")
    print(f"{dose_rate_Gy_s:.4g} Gy/s * {exposure_time_s} s = {cumulative_dose_Gy:.4g} Gy")
    
    return cumulative_dose_Gy

# Run the calculation and store the final answer
final_answer = calculate_cumulative_dose()

# The final answer is required in a specific format
print(f"\n<<< {final_answer:.4g} >>>")
