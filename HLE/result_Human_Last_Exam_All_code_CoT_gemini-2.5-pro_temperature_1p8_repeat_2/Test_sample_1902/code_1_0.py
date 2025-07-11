import math

def calculate_cumulative_dose():
    """
    Calculates the cumulative surface dose based on ionization chamber measurements.
    """
    # 1. Define given constants and physical parameters
    I = 2.0e-12         # Ionization chamber current in Amperes (C/s)
    beam_h_mm = 0.3       # Photon beam horizontal size at focus in mm
    beam_v_mm = 6.0       # Photon beam vertical size at focus in mm
    L_cm = 15.1         # Length of ionization chamber in cm
    rho_air_mg_cm3 = 1.293# Density of air in mg/cm^3
    t_exp_s = 0.02        # Effective exposure time for a point on the surface in seconds

    # Physical constant
    W_air_over_e = 33.97  # Average energy to create an ion pair in air, in Joules per Coulomb (J/C)

    # 2. Convert units to SI for consistency (meters, kilograms, seconds)
    beam_h_m = beam_h_mm / 1000.0
    beam_v_m = beam_v_mm / 1000.0
    L_m = L_cm / 100.0
    # 1 mg/cm^3 = 1 kg/m^3
    rho_air_kg_m3 = rho_air_mg_cm3

    # 3. Calculate the mass of the irradiated air in the chamber
    # Beam area in square meters
    A_m2 = beam_h_m * beam_v_m
    # Volume of irradiated air in cubic meters
    V_m3 = A_m2 * L_m
    # Mass of irradiated air in kilograms
    m_air_kg = V_m3 * rho_air_kg_m3

    # 4. Calculate the absorbed dose rate in Gray/s (J/kg/s)
    # Energy deposited per second (Power) in Joules/second
    P_dep_Js = I * W_air_over_e
    # Dose rate = Power / mass. Assumed D_rate_tissue = D_rate_air.
    dose_rate_Gy_s = P_dep_Js / m_air_kg

    # 5. Calculate the cumulative surface dose in Gray (Gy)
    cumulative_dose_Gy = dose_rate_Gy_s * t_exp_s

    # 6. Output the results
    print("The final equation for the cumulative dose is:")
    print("Dose = (Ionization Current * W_air/e * Exposure Time) / (Beam Area * Chamber Length * Air Density)")
    print("\nSubstituting the values (in SI units):")
    # Using more explicit variable names for clarity in the equation printout
    ionization_current_A = I
    beam_area_m2 = A_m2
    chamber_length_m = L_m
    air_density_kg_m3 = rho_air_kg_m3
    exposure_time_s = t_exp_s
    
    print(f"Dose = ({ionization_current_A:.1e} C/s * {W_air_over_e:.2f} J/C * {exposure_time_s:.2f} s) / (({beam_area_m2:.1e} m^2) * ({chamber_length_m:.3f} m) * ({air_density_kg_m3:.3f} kg/m^3))")

    # Convert final answer to milliGray (mGy) for better readability
    cumulative_dose_mGy = cumulative_dose_Gy * 1000
    
    # Format final answer to 3 significant figures
    final_answer_mGy = float(f"{cumulative_dose_mGy:.3g}")
    
    print(f"\nThe cumulative surface dose is {cumulative_dose_mGy:.6f} mGy.")
    print(f"Rounded to three significant figures, the dose is {final_answer_mGy} mGy.")
    
    # Final answer in the specified format
    print(f"<<<{final_answer_mGy}>>>")

if __name__ == "__main__":
    calculate_cumulative_dose()