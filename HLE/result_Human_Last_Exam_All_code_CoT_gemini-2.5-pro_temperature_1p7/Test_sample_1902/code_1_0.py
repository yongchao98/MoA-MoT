import math

def calculate_surface_dose():
    """
    Calculates the cumulative surface dose based on synchrotron experimental data.
    """

    # --- Input Parameters from the problem statement ---
    # Using SI units (meters, kilograms, seconds, Coulombs, etc.) for calculations
    beam_width_focus_mm = 0.3
    beam_height_mm = 6.0
    chamber_current_pa = 2.0
    air_density_mg_cm3 = 1.293
    chamber_length_cm = 15.1
    exposure_time_s = 0.02
    
    # --- Physical Constants ---
    W_air_ev = 33.97  # Mean energy to create an ion pair in air (eV/ip)
    e_coulomb = 1.602e-19  # Elementary charge in Coulombs (C/ip)
    eV_to_J = 1.602e-19  # Conversion factor from eV to Joules

    # --- Unit Conversions for Calculation ---
    beam_width_focus_m = beam_width_focus_mm / 1000.0
    beam_height_m = beam_height_mm / 1000.0
    chamber_current_a = chamber_current_pa * 1e-12
    # 1 mg/cm^3 = 1 (1e-6 kg) / (1e-6 m^3) = 1 kg/m^3
    air_density_kg_m3 = air_density_mg_cm3 
    chamber_length_m = chamber_length_cm / 100.0

    print("--- Cumulative Surface Dose Calculation ---")
    print("This script calculates the cumulative surface dose based on the provided experimental data.")
    print("\nMy Plan:")
    print("1. Calculate the mass of the air in the ionization chamber that is irradiated by the photon beam.")
    print("2. Calculate the energy deposited per second (Power) from the measured ionization current.")
    print("3. Calculate the absorbed dose rate (Power per unit mass).")
    print("4. Multiply the dose rate by the exposure time to find the cumulative dose.")
    
    # --- Step-by-Step Calculation ---

    # Step 1: Calculate the mass of the irradiated air (m_air)
    beam_area_m2 = beam_width_focus_m * beam_height_m
    irradiated_volume_m3 = beam_area_m2 * chamber_length_m
    irradiated_mass_kg = irradiated_volume_m3 * air_density_kg_m3

    # Step 2: Calculate the energy deposited per second (Power) in Watts (J/s)
    ion_pairs_per_sec = chamber_current_a / e_coulomb
    power_watts = ion_pairs_per_sec * W_air_ev * eV_to_J
    
    # Step 3: Calculate the dose rate in Gray/second (Gy/s)
    # Dose Rate = Energy per second / mass = (J/s) / kg = Gy/s
    dose_rate_gy_s = power_watts / irradiated_mass_kg

    # Step 4: Calculate the cumulative dose in Gray (Gy)
    cumulative_dose_gy = dose_rate_gy_s * exposure_time_s
    
    # --- Present the final equation with all numbers plugged in ---
    # The final equation combines all steps:
    # Dose = ( (Current / e) * W_air_in_Joules ) / ( Beam_Area * Chamber_Length * Air_Density ) * Exposure_Time
    print("\nThe final equation for cumulative dose is constructed by combining these steps.")
    print("Plugging in the numerical values:\n")
    print(f"Cumulative Dose (Gy) = ( (({chamber_current_pa}e-12 A) / ({e_coulomb} C/ip)) * ({W_air_ev} eV/ip) * ({eV_to_J} J/eV) ) / "
          f"( (({beam_width_focus_mm}/1000) m * ({beam_height_mm}/1000) m) * ({chamber_length_cm}/100) m * ({air_density_kg_m3}) kg/m^3 ) * ({exposure_time_s}) s")
    
    # --- Final Result ---
    # Convert Gy to mGy for a more readable output
    cumulative_dose_mgy = cumulative_dose_gy * 1000
    
    print("\n--- Result ---")
    print(f"The calculated cumulative surface dose is: {cumulative_dose_gy:.4e} Gy")
    print(f"This is equivalent to {cumulative_dose_mgy:.5f} mGy.")

    # Return final answer in the required format
    print(f"\n<<<{cumulative_dose_mgy:.5f}>>>")

calculate_surface_dose()