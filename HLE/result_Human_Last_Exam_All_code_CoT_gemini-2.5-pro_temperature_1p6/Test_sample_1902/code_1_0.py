import math

def calculate_surface_dose():
    """
    Calculates the cumulative surface dose based on synchrotron imaging parameters.
    """
    # --- Given Parameters ---
    beam_width_focus_mm = 0.3  # mm
    beam_height_focus_mm = 6.0   # mm
    ionization_current_pA = 2.0  # pA
    air_density_mg_cm3 = 1.293  # mg/cm^3
    chamber_length_cm = 15.1  # cm
    mass_energy_absorption_coeff_air_cm2_g = 0.328  # cm^2/g
    total_exposure_time_s = 0.02  # s

    # --- Physical Constants ---
    # W_air/e: Average energy to create an ion pair in air / elementary charge
    W_over_e_air_J_C = 33.97  # J/C

    # --- Unit Conversions for a consistent system (cm, g, s, C) ---
    beam_width_focus_cm = beam_width_focus_mm / 10.0
    beam_height_focus_cm = beam_height_focus_mm / 10.0
    ionization_current_A = ionization_current_pA * 1e-12  # A is C/s
    air_density_g_cm3 = air_density_mg_cm3 / 1000.0  # g/cm^3

    # --- Calculations ---
    
    # Calculate beam area at focus in cm^2
    beam_area_cm2 = beam_width_focus_cm * beam_height_focus_cm

    # Calculate total power absorbed in the ion chamber (Joules/sec)
    power_absorbed_J_s = ionization_current_A * W_over_e_air_J_C

    # Calculate the exponent for the absorption formula: mu_en * L
    exponent = mass_energy_absorption_coeff_air_cm2_g * air_density_g_cm3 * chamber_length_cm

    # Calculate the fraction of the beam's energy absorbed by the air in the chamber
    fraction_absorbed = 1 - math.exp(-exponent)

    # Step 1: Calculate Energy Fluence Rate (Psi_dot) in J/(s*cm^2)
    psi_dot_J_s_cm2 = power_absorbed_J_s / (beam_area_cm2 * fraction_absorbed)

    # Step 2: Calculate Dose Rate (D_dot)
    # Dose rate in J/(s*g)
    dose_rate_J_s_g = psi_dot_J_s_cm2 * mass_energy_absorption_coeff_air_cm2_g
    # Convert dose rate to Gy/s (which is J/(s*kg))
    dose_rate_Gy_s = dose_rate_J_s_g * 1000.0

    # Step 3: Calculate Cumulative Dose (D_cumulative)
    cumulative_dose_Gy = dose_rate_Gy_s * total_exposure_time_s

    # --- Output the step-by-step calculation ---
    print("Calculation of the cumulative surface dose:\n")

    print("Step 1: Calculate the energy fluence rate (Ψ_dot) at the focus.")
    print("Formula: Ψ_dot = (I * (W/e)) / (A * (1 - exp(-(μ_en/ρ) * ρ * L)))")
    print(f"Ψ_dot = ({ionization_current_A:.1e} C/s * {W_over_e_air_J_C} J/C) / (({beam_width_focus_cm} cm * {beam_height_focus_cm} cm) * (1 - exp(-{mass_energy_absorption_coeff_air_cm2_g} cm^2/g * {air_density_g_cm3:.6f} g/cm^3 * {chamber_length_cm} cm)))")
    print(f"Ψ_dot = {power_absorbed_J_s:.3e} J/s / ({beam_area_cm2:.3f} cm^2 * {fraction_absorbed:.5f})")
    print(f"Result: Ψ_dot = {psi_dot_J_s_cm2:.3e} J/s/cm^2\n")

    print("Step 2: Calculate the dose rate (D_dot) in tissue.")
    print("Assumption: Tissue absorption is approximately equal to air absorption.")
    print("Formula: D_dot = Ψ_dot * (μ_en/ρ)")
    print(f"D_dot = {psi_dot_J_s_cm2:.3e} J/s/cm^2 * {mass_energy_absorption_coeff_air_cm2_g} cm^2/g * 1000 g/kg")
    print(f"Result: D_dot = {dose_rate_Gy_s:.5f} Gy/s\n")

    print("Step 3: Calculate the cumulative surface dose (D_cumulative).")
    print("Formula: D_cumulative = D_dot * t_exposure")
    print(f"D_cumulative = {dose_rate_Gy_s:.5f} Gy/s * {total_exposure_time_s} s")
    print(f"Result: D_cumulative = {cumulative_dose_Gy:.3e} Gy")
    print(f"\nThe cumulative surface dose is approximately {cumulative_dose_Gy * 1e6:.2f} microGray (µGy).")
    
    # Return the final numeric value for the answer tag
    return cumulative_dose_Gy

# Execute the calculation and store the result
final_dose = calculate_surface_dose()
print(f"\n<<<3.88e-6>>>")