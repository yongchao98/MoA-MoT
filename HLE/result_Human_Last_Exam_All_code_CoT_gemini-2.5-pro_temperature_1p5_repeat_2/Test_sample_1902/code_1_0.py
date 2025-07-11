import math

def calculate_surface_dose():
    """
    Calculates the cumulative surface dose to a tissue based on synchrotron imaging parameters.
    """
    # --- Given Information & Constants ---
    # Using SI units for all calculations (meters, kg, seconds, Joules, Coulombs)

    # Beam and chamber dimensions
    beam_h_focus_mm = 0.3  # mm
    beam_v_focus_mm = 6.0  # mm
    chamber_length_cm = 15.1 # cm

    # Beam and material properties
    photon_energy_keV = 33.0  # keV
    chamber_current_pA = 2.0  # pA
    air_density_mg_cm3 = 1.293 # mg/cm^3
    mu_en_rho_air_cm2_g = 0.328 # cm^2/g (mass energy absorption coefficient of air)

    # Exposure time parameter
    exposure_time_s = 0.02 # s (dwell time for a point on the sample)

    # Physical constants
    W_air_eV_per_ip = 33.97 # eV/ion pair (average energy to create an ion pair in air)
    # The value W_air / e_charge is 33.97 J/C
    e_charge_C = 1.60217663e-19 # Coulombs (elementary charge)

    # --- Unit Conversions to SI ---
    beam_h_focus_m = beam_h_focus_mm / 1000.0
    beam_v_focus_m = beam_v_focus_mm / 1000.0
    chamber_length_m = chamber_length_cm / 100.0
    chamber_current_A = chamber_current_pA * 1e-12
    # 1 mg/cm^3 = 1 kg/m^3
    air_density_kg_m3 = air_density_mg_cm3 * 1.0
    # 1 cm^2/g = 0.1 m^2/kg
    mu_en_rho_air_m2_kg = mu_en_rho_air_cm2_g * 0.1
    photon_energy_J = photon_energy_keV * 1000 * e_charge_C
    W_air_J_per_C = W_air_eV_per_ip

    # --- Step-by-step Calculation ---

    # 1. Calculate the energy deposited per second (Power) in the ionization chamber (J/s)
    power_deposited_chamber_J_s = chamber_current_A * W_air_J_per_C

    # 2. Calculate the fraction of the beam's energy deposited in the chamber
    # First, find the linear energy absorption coefficient of air (mu_en)
    mu_en_air_per_m = mu_en_rho_air_m2_kg * air_density_kg_m3
    # Then, find the dimensionless exponent for the absorption formula
    exponent = mu_en_air_per_m * chamber_length_m
    # Finally, calculate the fraction of energy deposited
    fraction_deposited = 1 - math.exp(-exponent)

    # 3. Calculate the total photon flux (photons/s) entering the chamber
    total_energy_flux_J_s = power_deposited_chamber_J_s / fraction_deposited
    total_photon_flux_photons_s = total_energy_flux_J_s / photon_energy_J

    # 4. Calculate the photon fluence rate at the focus (photons/s/m^2)
    beam_area_focus_m2 = beam_h_focus_m * beam_v_focus_m
    fluence_rate_photons_s_m2 = total_photon_flux_photons_s / beam_area_focus_m2

    # 5. Calculate the dose rate at the surface of the tissue (Gy/s)
    # Assume (mu_en/rho)_tissue is approx. (mu_en/rho)_air
    mu_en_rho_tissue_m2_kg = mu_en_rho_air_m2_kg
    dose_rate_Gy_s = fluence_rate_photons_s_m2 * photon_energy_J * mu_en_rho_tissue_m2_kg

    # 6. Calculate the cumulative surface dose (Gy)
    cumulative_dose_Gy = dose_rate_Gy_s * exposure_time_s
    cumulative_dose_mGy = cumulative_dose_Gy * 1000

    # --- Output the results ---
    print("--- Calculation of Cumulative Surface Dose ---")

    print("\nStep 1: Calculate power deposited in chamber")
    print(f"Power = Chamber Current (A) * (W/e) for air (J/C)")
    print(f"Power = {chamber_current_A:.2e} A * {W_air_J_per_C:.2f} J/C = {power_deposited_chamber_J_s:.3e} J/s")

    print("\nStep 2: Calculate photon flux entering the chamber")
    print(f"Fraction of energy deposited = 1 - exp(-({mu_en_air_per_m:.4f} m^-1 * {chamber_length_m:.3f} m)) = {fraction_deposited:.5f}")
    print(f"Total photon flux = ({power_deposited_chamber_J_s:.3e} J/s) / ({fraction_deposited:.5f} * {photon_energy_J:.3e} J/photon) = {total_photon_flux_photons_s:.3e} photons/s")

    print("\nStep 3: Calculate dose rate at tissue surface")
    print(f"Photon fluence rate = {total_photon_flux_photons_s:.3e} photons/s / {beam_area_focus_m2:.2e} m^2 = {fluence_rate_photons_s_m2:.3e} photons/s/m^2")
    print(f"Dose Rate = Fluence Rate * Photon Energy (J) * Mass Energy Absorption Coeff. (m^2/kg)")
    print(f"Dose Rate = {fluence_rate_photons_s_m2:.3e} photons/s/m^2 * {photon_energy_J:.3e} J * {mu_en_rho_tissue_m2_kg:.4f} m^2/kg = {dose_rate_Gy_s:.4f} Gy/s")

    print("\nStep 4: Calculate final cumulative dose")
    print(f"Cumulative Dose (Gy) = Dose Rate (Gy/s) * Exposure Time (s)")
    print(f"Cumulative Dose = {dose_rate_Gy_s:.4f} Gy/s * {exposure_time_s:.2f} s = {cumulative_dose_Gy:.6f} Gy")

    print(f"\nFinal Cumulative Surface Dose = {cumulative_dose_mGy:.2f} mGy")
    
    # Return the value for the final answer tag
    return cumulative_dose_mGy

if __name__ == '__main__':
    final_answer = calculate_surface_dose()
    # The final answer is printed within the function.
    # To conform to the format, the raw number is needed for the tag.
    # print(f"\n<<< {final_answer:.2f} >>>") # For generating the tag.

calculate_surface_dose()