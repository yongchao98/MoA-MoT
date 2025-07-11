import math

def calculate_surface_dose():
    """
    Calculates the cumulative surface dose based on synchrotron experimental parameters.
    """
    # Step 1: Define constants and given values from the problem description.
    # Beam parameters
    beam_width_h_mm = 0.3  # mm
    beam_width_v_mm = 6.0    # mm

    # Ionization chamber parameters
    current_pA = 2.0  # pA
    chamber_length_cm = 15.1  # cm

    # Physical constants and material properties
    air_density_mg_cm3 = 1.293  # mg/cm^3
    W_air_over_e_J_C = 33.97  # J/C, mean energy to create an ion pair in dry air

    # Scan parameter (interpreted as the total exposure time for a point on the surface)
    total_exposure_time_s = 0.02  # s

    # Step 2: Convert units to a consistent system (cm, g, s and then kg for dose).
    # Convert beam dimensions from mm to cm
    beam_width_h_cm = beam_width_h_mm / 10.0
    beam_width_v_cm = beam_width_v_mm / 10.0

    # Convert current from pA to A (C/s)
    current_A = current_pA * 1e-12

    # Convert air density from mg/cm^3 to g/cm^3
    air_density_g_cm3 = air_density_mg_cm3 / 1000.0

    # Step 3: Calculate the mass of the irradiated air.
    # Calculate beam area at focus in cm^2
    beam_area_cm2 = beam_width_h_cm * beam_width_v_cm

    # Calculate irradiated volume of the ionization chamber in cm^3
    irradiated_volume_cm3 = beam_area_cm2 * chamber_length_cm

    # Calculate the mass of the irradiated air in grams, then kilograms.
    air_mass_g = irradiated_volume_cm3 * air_density_g_cm3
    air_mass_kg = air_mass_g / 1000.0

    # Step 4: Calculate the dose rate in Gy/s (J/kg/s).
    # The problem assumes tissue absorption properties are similar to air.
    # Energy deposited per second (J/s) = Current (C/s) * Energy per charge (J/C)
    energy_deposited_per_second_J_s = current_A * W_air_over_e_J_C
    
    # Dose rate (Gy/s) = Energy deposited per second (J/s) / mass (kg)
    dose_rate_Gy_s = energy_deposited_per_second_J_s / air_mass_kg

    # Step 5: Calculate the total cumulative dose.
    # Total Dose (Gy) = Dose Rate (Gy/s) * Total Exposure Time (s)
    total_dose_Gy = dose_rate_Gy_s * total_exposure_time_s

    # Convert total dose to milliGrays (mGy) for a more convenient representation.
    total_dose_mGy = total_dose_Gy * 1000

    # Step 6: Print the final equation and the result.
    print("The final calculation for the cumulative surface dose is structured as follows:")
    print("Dose = (Dose Rate) * (Exposure Time)")
    print("where Dose Rate = (Current * W_air/e) / (Beam Area * Chamber Length * Air Density)")

    print("\n--- Final Equation with a-l Numbers ---")
    # This equation shows all the initial values and conversion factors used.
    equation_str = (
        f"Dose (mGy) = (({current_pA}e-12 A * {W_air_over_e_J_C} J/C) / "
        f"(({beam_width_h_mm}/10 * {beam_width_v_mm}/10) cm^2 * {chamber_length_cm} cm * {air_density_mg_cm3}/1000 g/cm^3 * 0.001 kg/g))"
        f" * {total_exposure_time_s} s * 1000 mGy/Gy"
    )
    print(equation_str)

    # Print the final numerical answer
    print(f"\nThe cumulative surface dose is: {total_dose_mGy:.5f} mGy")


# Execute the calculation
calculate_surface_dose()