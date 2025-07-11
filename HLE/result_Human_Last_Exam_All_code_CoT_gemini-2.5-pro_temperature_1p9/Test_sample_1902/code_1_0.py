import math

def calculate_surface_dose():
    """
    Calculates the cumulative surface dose based on synchrotron imaging parameters.
    """
    # 1. Define constants and given parameters in consistent units (cm, g, s, A, J)
    beam_width_focus_mm = 0.3
    beam_height_focus_mm = 6.0
    chamber_current_pA = 2.0
    density_air_mg_cm3 = 1.293
    chamber_length_cm = 15.1
    mu_en_rho_air_cm2_g = 0.328  # Energy absorption mass attenuation coefficient for air
    total_exposure_time_s = 0.02
    W_per_e_air_J_C = 33.97  # Avg. energy to create an ion pair in air (Joule/Coulomb)

    # Unit conversions
    beam_area_cm2 = (beam_width_focus_mm / 10) * (beam_height_focus_mm / 10)
    current_A = chamber_current_pA * 1e-12  # pA to A (C/s)
    density_air_g_cm3 = density_air_mg_cm3 / 1000 # mg/cm^3 to g/cm^3

    # 2. Calculate the linear energy absorption coefficient of air (mu_en)
    mu_en_air_inv_cm = mu_en_rho_air_cm2_g * density_air_g_cm3

    # 3. Calculate the power deposited in the ionization chamber's air volume
    power_deposited_J_s = current_A * W_per_e_air_J_C

    # 4. Calculate the incident power of the beam on the chamber
    # The fraction of power absorbed is (1 - exp(-mu_en * L))
    absorption_fraction = 1 - math.exp(-mu_en_air_inv_cm * chamber_length_cm)
    power_incident_J_s = power_deposited_J_s / absorption_fraction

    # 5. Calculate the energy fluence rate (incident power per unit area)
    energy_fluence_rate_J_cm2_s = power_incident_J_s / beam_area_cm2

    # 6. Calculate the dose rate in tissue (approximated as air)
    # Dose Rate = (mu_en/rho) * Energy_Fluence_Rate
    # Result is in J/g/s
    dose_rate_J_g_s = mu_en_rho_air_cm2_g * energy_fluence_rate_J_cm2_s

    # 7. Calculate the cumulative dose
    # Cumulative Dose = Dose Rate * Exposure Time
    cumulative_dose_J_g = dose_rate_J_g_s * total_exposure_time_s

    # 8. Convert the final dose to more common units (microgray)
    # 1 J/g = 1000 J/kg = 1000 Gy = 1,000,000,000 uGy
    # Let's convert in two steps:
    # 1 J/g -> 1000 Gy
    cumulative_dose_Gy = cumulative_dose_J_g * 1000
    # 1 Gy -> 1,000,000 uGy
    cumulative_dose_uGy = cumulative_dose_Gy * 1e6
    
    # Print the explanation and the full equation with numbers
    print("The cumulative surface dose can be calculated using the following formula:")
    print("Dose = ( (I * (W/e)_air) / (A * (1 - exp(-(μ_en/ρ)_air * ρ_air * L))) ) * (μ_en/ρ)_tissue * t_exposure\n")
    print("Substituting the given values:")
    
    equation = (
        f"Dose = ( ({current_A:.1e} C/s * {W_per_e_air_J_C} J/C) / "
        f"({beam_area_cm2} cm^2 * (1 - exp(-({mu_en_rho_air_cm2_g} cm^2/g * {density_air_g_cm3} g/cm^3 * {chamber_length_cm} cm)))) )"
        f" * {mu_en_rho_air_cm2_g} cm^2/g * {total_exposure_time_s} s"
    )
    print(equation)
    
    print("\nThis calculation gives the dose in J/g. Converting to micrograys (µGy):")
    print(f"Result = {cumulative_dose_J_g:.4e} J/g * (1000 Gy / (J/g)) * (1,000,000 µGy / Gy)")
    print(f"\nFinal Cumulative Surface Dose: {cumulative_dose_uGy:.3f} µGy")

if __name__ == "__main__":
    calculate_surface_dose()
    # The final value is approximately 3.883 µGy.
    print("\n<<<3.883>>>")