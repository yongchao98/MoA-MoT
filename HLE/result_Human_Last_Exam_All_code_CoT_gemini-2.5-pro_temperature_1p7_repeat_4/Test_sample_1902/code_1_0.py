import math

def calculate_surface_dose():
    """
    Calculates the cumulative surface dose based on synchrotron imaging system parameters.
    """
    # --- Input Parameters ---
    # Beam and Chamber Geometry
    beam_width_h_focus_mm = 0.3  # mm
    beam_height_v_mm = 6.0       # mm
    chamber_length_cm = 15.1     # cm

    # Electrical and Physical Properties
    current_pA = 2.0             # pA
    rho_air_mg_cm3 = 1.293       # mg/cm^3
    mu_en_rho_air_cm2_g = 0.328  # cm^2/g

    # Time and Energy Constants
    dwell_time_s = 0.02          # s
    # Average energy to create an ion pair in air (W_air)
    W_air_J_per_C = 33.97        # J/C

    # --- Unit Conversions for Calculation ---
    # Convert all inputs to a consistent unit system (cm, g, s, A, C)
    beam_width_h_focus_cm = beam_width_h_focus_mm / 10.0
    beam_height_v_cm = beam_height_v_mm / 10.0
    current_A = current_pA * 1e-12
    rho_air_g_cm3 = rho_air_mg_cm3 / 1000.0

    # --- Step 1: Calculate the Dose Rate (dD/dt) ---
    # The dose rate (dD/dt) at the surface is given by Ψ * (μ_en/ρ)_tissue
    # where Ψ is the energy fluence rate. We can derive a formula relating this to the
    # ionization chamber current.
    # dD/dt [J/(g*s)] = (I * W_air * (μ_en/ρ)_air) / (A_beam * (1 - exp(-(μ_en/ρ)_air * ρ_air * L)))
    
    # Calculate beam area at focus
    A_beam_cm2 = beam_width_h_focus_cm * beam_height_v_cm

    # Calculate the exponent term for the absorption calculation
    # This term must be dimensionless: (cm^2/g) * (g/cm^3) * cm
    exponent = mu_en_rho_air_cm2_g * rho_air_g_cm3 * chamber_length_cm

    # Calculate the dose rate in J/(g*s)
    numerator = current_A * W_air_J_per_C * mu_en_rho_air_cm2_g
    denominator = A_beam_cm2 * (1 - math.exp(-exponent))
    dose_rate_J_g_s = numerator / denominator

    # Convert dose rate from J/g*s to Gy/s (1 J/kg = 1 Gy, so 1 J/g = 1000 Gy)
    dose_rate_Gy_s = dose_rate_J_g_s * 1000.0

    # --- Step 2: Calculate Cumulative Dose ---
    # Cumulative Dose = Dose Rate * Dwell Time
    cumulative_dose_Gy = dose_rate_Gy_s * dwell_time_s

    # --- Print the Calculation Breakdown and Final Answer ---
    print("--- Dose Rate Calculation ---")
    print("The equation for Dose Rate is:")
    print("Dose_Rate [Gy/s] = (I * W_air * (μ_en/ρ)_air) / (A_beam * (1 - exp(-(μ_en/ρ)_air * ρ_air * L))) * 1000\n")
    print("Plugging in the numbers:")
    print(f"Dose_Rate [Gy/s] = ({current_A:.2e} A * {W_air_J_per_C:.2f} J/C * {mu_en_rho_air_cm2_g:.3f} cm^2/g) / "
          f"({A_beam_cm2:.3f} cm^2 * (1 - exp(-{mu_en_rho_air_cm2_g:.3f} * {rho_air_g_cm3:.5f} * {chamber_length_cm:.1f}))) * 1000")
    print(f"Dose_Rate = {dose_rate_Gy_s:.6f} Gy/s\n")
    
    print("--- Cumulative Dose Calculation ---")
    print("The equation for Cumulative Dose is:")
    print("Cumulative_Dose [Gy] = Dose_Rate [Gy/s] * Dwell_Time [s]\n")
    print("Plugging in the numbers:")
    print(f"Cumulative_Dose [Gy] = {dose_rate_Gy_s:.6f} Gy/s * {dwell_time_s:.2f} s")
    
    # Final result
    print(f"\nFinal Calculated Cumulative Surface Dose = {cumulative_dose_Gy:.3e} Gy")
    
    return cumulative_dose_Gy

# Run the calculation
final_dose = calculate_surface_dose()
# The final answer block is below this script, formatted as requested.
# print(f"<<<{final_dose:.3e}>>>") # This is for internal check

calculate_surface_dose()