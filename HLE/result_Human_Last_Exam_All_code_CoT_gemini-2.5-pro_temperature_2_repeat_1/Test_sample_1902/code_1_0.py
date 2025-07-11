import math

def calculate_cumulative_dose():
    """
    Calculates the cumulative surface dose based on synchrotron experimental parameters.
    """
    # Step 1: Define the given parameters and constants in SI units.
    
    # Ionization chamber current in Amperes (C/s)
    I = 2.0e-12  # A (from 2.0 pA)
    
    # Photon beam dimensions at focus in meters
    beam_width_focus = 0.3e-3  # m (from 0.3 mm)
    beam_height = 6.0e-3     # m (from 6.0 mm)
    
    # Length of the ionization chamber in meters
    L_chamber = 15.1e-2  # m (from 15.1 cm)
    
    # Density of air in kg/m^3
    # Given: 1.293 mg/cm^3 = 1.293 * (1e-6 kg) / (1e-2 m)^3 = 1.293 kg/m^3
    rho_air = 1.293  # kg/m^3
    
    # Average energy to create an ion pair in air, in Joules per Coulomb (J/C)
    W_air_per_charge = 33.97  # J/C
    
    # Exposure time for any point on the surface, in seconds
    exposure_time = 0.02  # s

    # Step 2: Calculate the rate of energy deposited (Power) in the chamber air.
    # Power (J/s) = Current (C/s) * Energy per charge (J/C)
    power_deposited = I * W_air_per_charge
    
    # Step 3: Calculate the mass of the air being irradiated.
    # Volume (m^3) = width (m) * height (m) * length (m)
    beam_area = beam_width_focus * beam_height
    air_volume = beam_area * L_chamber
    # Mass (kg) = Volume (m^3) * Density (kg/m^3)
    air_mass = air_volume * rho_air

    # Step 4: Calculate the dose rate in Gray/second (J/kg/s).
    # We assume the dose rate in tissue is equal to the dose rate in air.
    # Dose Rate (Gy/s) = Power (J/s) / Mass (kg)
    dose_rate = power_deposited / air_mass
    
    # Step 5: Calculate the final cumulative dose in Gray (Gy).
    # Cumulative Dose (Gy) = Dose Rate (Gy/s) * Exposure Time (s)
    cumulative_dose = dose_rate * exposure_time
    
    # Step 6: Print the equations and results as requested.
    print("This script calculates the cumulative surface dose to a tissue.")
    print("-" * 60)
    
    print("Calculation of the dose rate:")
    print(f"Dose Rate = (Deposited Power) / (Mass of Air)")
    print(f"Dose Rate = ({power_deposited:.3e} J/s) / ({air_mass:.3e} kg)")
    print(f"Dose Rate = {dose_rate:.3e} Gy/s\n")
    
    print("Calculation of the cumulative dose:")
    print(f"Cumulative Dose = Dose Rate * Exposure Time")
    print(f"Cumulative Dose = {dose_rate:.3e} Gy/s * {exposure_time} s")
    print(f"Cumulative Dose = {cumulative_dose:.3e} Gy")
    
    print("-" * 60)
    print(f"The final calculated cumulative surface dose is {cumulative_dose} Gy.")


calculate_cumulative_dose()