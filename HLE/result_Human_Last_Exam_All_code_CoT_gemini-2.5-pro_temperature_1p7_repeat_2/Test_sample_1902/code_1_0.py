import math

def calculate_surface_dose():
    """
    Calculates the cumulative surface dose to a tissue based on ionization chamber measurements.
    """
    # --- Given Information ---
    # Photon beam size at focus (horizontal by vertical)
    beam_width_h_mm = 0.3  # mm
    beam_width_v_mm = 6.0    # mm

    # Ionization chamber parameters
    I_pA = 2.0  # pA
    rho_air_mg_cm3 = 1.293  # mg/cm^3
    L_cm = 15.1  # cm

    # Other parameters
    t_eff_s = 0.02  # s (Effective exposure time for a point on the surface)

    # --- Physical Constants ---
    # Mean energy to create an ion pair in dry air, per unit charge (J/C)
    W_air_over_e = 33.97

    # --- Unit Conversions ---
    # Convert beam dimensions from mm to cm
    beam_width_h_cm = beam_width_h_mm / 10
    beam_width_v_cm = beam_width_v_mm / 10

    # Convert current from pA to A (C/s)
    I_A = I_pA * 1e-12

    # Convert air density from mg/cm^3 to g/cm^3 and kg/m^3
    rho_air_g_cm3 = rho_air_mg_cm3 / 1000
    rho_air_kg_cm3 = rho_air_g_cm3 / 1000

    # --- Step 1: Calculate the irradiated volume and mass of air ---
    # Beam area in cm^2
    beam_area_cm2 = beam_width_h_cm * beam_width_v_cm

    # Irradiated volume in the chamber in cm^3
    volume_cm3 = beam_area_cm2 * L_cm

    # Mass of the irradiated air in kg
    mass_air_kg = volume_cm3 * rho_air_kg_cm3
    
    # --- Step 2: Calculate the Dose Rate ---
    # Dose Rate in Gy/s (J/kg/s)
    # Formula: Dose Rate = (Current / mass) * (W_air / e)
    dose_rate_Gy_s = (I_A / mass_air_kg) * W_air_over_e
    
    # --- Step 3: Calculate the Cumulative Dose ---
    # Cumulative Dose in Gy
    # Formula: Cumulative Dose = Dose Rate * Effective Exposure Time
    cumulative_dose_Gy = dose_rate_Gy_s * t_eff_s

    # --- Print the output as requested ---
    print("This script calculates the cumulative surface dose to the tissue.")
    print("-" * 60)
    print("Calculation Steps:")
    
    # Print mass calculation
    print("1. Calculate the mass of the irradiated air (m).")
    print(f"   m = (Beam Area at Focus * Chamber Length) * Air Density")
    print(f"   m = (({beam_width_h_cm} cm * {beam_width_v_cm} cm) * {L_cm} cm) * {rho_air_kg_cm3:.4e} kg/cm^3")
    print(f"   m = {mass_air_kg:.4e} kg")
    print()

    # Print dose rate calculation
    print("2. Calculate the dose rate to air (D_dot).")
    print(f"   D_dot = (Ionization Current / Mass of Air) * (W_air / e)")
    print(f"   D_dot = ({I_A:.1e} C/s / {mass_air_kg:.4e} kg) * {W_air_over_e} J/C")
    print(f"   D_dot = {dose_rate_Gy_s:.4e} Gy/s")
    print()

    # Print cumulative dose calculation
    print("3. Calculate the cumulative surface dose (D).")
    print(f"   D = Dose Rate * Effective Exposure Time")
    print(f"   D = {dose_rate_Gy_s:.4e} Gy/s * {t_eff_s} s")
    print(f"   D = {cumulative_dose_Gy:.4e} Gy")
    print("-" * 60)
    
    # Return the final numerical answer for the submission format
    return cumulative_dose_Gy

# Run the calculation and store the result
final_dose = calculate_surface_dose()

# The final answer in the required format
print(f"<<<{final_dose:.3e}>>>")
