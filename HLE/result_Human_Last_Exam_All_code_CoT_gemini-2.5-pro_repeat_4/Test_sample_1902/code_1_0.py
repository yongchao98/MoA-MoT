import math

def calculate_surface_dose():
    """
    Calculates the cumulative surface dose based on ionization chamber measurements.
    """
    # --- Given Information (with unit conversions to SI) ---

    # Ionization chamber current in Amperes (C/s)
    # 2.0 pA = 2.0e-12 A
    current_I = 2.0e-12  # C/s

    # Photon beam size at focus (horizontal x vertical) in meters
    # 0.3 mm x 6 mm = 0.0003 m x 0.006 m
    beam_width_h_focus = 0.3 / 1000  # m
    beam_width_v_focus = 6 / 1000   # m
    beam_area_focus = beam_width_h_focus * beam_width_v_focus # m^2

    # Length of ionization chamber in meters
    # 15.1 cm = 0.151 m
    chamber_length_L = 15.1 / 100  # m

    # Density of air in kg/m^3
    # 1.293 mg/cm^3 = 1.293 * (1e-6 kg) / (1e-2 m)^3 = 1.293 kg/m^3
    density_air_rho = 1.293  # kg/m^3

    # Exposure time for a point on the surface in seconds
    # Given by the "Ratio of beam’s width at focus to the scan step of subject per exposure time"
    exposure_time_t = 0.02  # s

    # --- Physical Constants ---

    # Average energy to create an ion pair in air in Joules/Coulomb
    W_air = 33.97  # J/C

    # --- Calculation Steps ---

    # 1. Calculate the energy absorbed per second (Power) in the chamber
    power_absorbed = current_I * W_air

    # 2. Calculate the volume of air irradiated by the beam
    irradiated_volume = beam_area_focus * chamber_length_L

    # 3. Calculate the mass of the irradiated air
    mass_air = irradiated_volume * density_air_rho

    # 4. Calculate the dose rate in air (in Gray/second)
    # Dose is Energy/Mass (J/kg), so Dose Rate is Power/Mass (J/s / kg = Gy/s)
    dose_rate = power_absorbed / mass_air

    # 5. Calculate the cumulative dose
    # This is the dose rate multiplied by the exposure time.
    # We assume the dose to tissue is the same as the dose to air.
    cumulative_dose = dose_rate * exposure_time_t
    
    # --- Output the results ---
    
    print("This script calculates the cumulative surface dose to a tissue.")
    print("The calculation is based on the formula:")
    print("Cumulative Dose = ( (Current * W_air) / (Beam Area * Chamber Length * Air Density) ) * Exposure Time\n")
    
    print("Substituting the given values:")
    
    # Using a formatted string to show the equation with numbers
    # Each value is formatted using scientific notation for clarity
    equation_str = (
        f"Cumulative Dose = ( ({current_I:.1e} A * {W_air:.2f} J/C) / "
        f"({beam_area_focus:.1e} m^2 * {chamber_length_L:.3f} m * {density_air_rho:.3f} kg/m^3) ) * "
        f"{exposure_time_t:.2f} s"
    )
    print(equation_str)
    
    print(f"\nStep 1: Absorbed Power = {power_absorbed:.3e} J/s")
    print(f"Step 2: Irradiated Mass of Air = {mass_air:.3e} kg")
    print(f"Step 3: Dose Rate = {dose_rate:.3e} Gy/s")
    print(f"Step 4: Cumulative Dose = {cumulative_dose:.3e} Gy")
    
    # Returning the final answer for the wrapper
    return cumulative_dose

# Execute the calculation and print the final answer in the required format
final_dose = calculate_surface_dose()
print(f"\nThe calculated cumulative surface dose is {final_dose:.3g} Gy.")
print(f"<<<{final_dose:.3g}>>>")
