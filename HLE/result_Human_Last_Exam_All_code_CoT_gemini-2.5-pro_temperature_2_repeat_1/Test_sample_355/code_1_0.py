import numpy as np

def calculate_exposure_time():
    """
    Calculates the exposure time needed to measure a star's B-band magnitude
    to a given accuracy.
    """
    # 1. Calculate required Signal-to-Noise Ratio (S/N)
    delta_m = 0.01  # Desired magnitude uncertainty
    
    # The relation is delta_m ≈ (2.5 / ln(10)) * (1 / S/N)
    snr_factor = 2.5 / np.log(10)
    required_snr = snr_factor / delta_m
    
    print("Step 1: Calculate Required Signal-to-Noise Ratio (S/N)")
    print(f"The required S/N = (2.5 / ln(10)) / δm = {snr_factor:.4f} / {delta_m} = {required_snr:.2f}")
    print("-" * 50)

    # 2. Calculate total photons needed (N)
    # For Poisson statistics, S/N = sqrt(N) => N = (S/N)^2
    total_photons_needed = required_snr**2
    
    print("Step 2: Calculate Total Photons Needed (N)")
    print(f"Total photons needed (N) = S/N² = {required_snr:.2f}² = {total_photons_needed:.1f}")
    print("-" * 50)
    
    # 3. Calculate the rate of incoming photons (R)
    
    # Astrophysical and physical constants
    star_mag = 20.0
    telescope_diameter_m = 1.0
    # B-band zero-point flux (Bessell 2005) in cgs units
    f_lambda_0_cgs = 4.26e-9  # erg cm⁻² s⁻¹ Å⁻¹
    # B-band parameters
    lambda_b_nm = 440.0 # Central wavelength in nm
    delta_lambda_b_nm = 98.0 # Bandwidth in nm
    
    # Constants for conversion
    h_joule_s = 6.626e-34 # Planck's constant
    c_m_s = 3.0e8       # Speed of light
    
    # Convert B-band parameters to meters
    lambda_b_m = lambda_b_nm * 1e-9
    delta_lambda_b_m = delta_lambda_b_nm * 1e-9

    print("Step 3: Calculate the Rate of Incoming Photons (R)")
    
    # a) Calculate photon flux for a magnitude 0 star
    
    # Convert F_lambda from cgs (erg cm⁻² s⁻¹ Å⁻¹) to SI (W m⁻³ or J s⁻¹ m⁻³)
    # 1 erg = 1e-7 J; 1 cm^2 = 1e-4 m^2; 1 Å = 1e-10 m
    # Factor = (1e-7 J/erg) / (1e-4 m²/cm²) / (1e-10 m/Å) = 1e7
    f_lambda_0_si = f_lambda_0_cgs * 1e7
    
    # Calculate energy of a single B-band photon
    e_photon_j = (h_joule_s * c_m_s) / lambda_b_m
    
    # Calculate photon flux rate per area per bandwidth
    photon_flux_density_0 = f_lambda_0_si / e_photon_j
    
    # Integrate over the bandwidth to get photon rate per area
    photon_rate_per_area_0 = photon_flux_density_0 * delta_lambda_b_m
    
    print("  a) For a magnitude 0 star:")
    print(f"    Photon energy = (h * c) / λ = ({h_joule_s:.3e} * {c_m_s:.1e}) / {lambda_b_m:.3e} = {e_photon_j:.3e} J")
    print(f"    Photon flux rate (mag 0) = F₀ / E_photon * Δλ = {f_lambda_0_cgs:.2e} * 1e7 / {e_photon_j:.3e} * {delta_lambda_b_m:.3e} = {photon_rate_per_area_0:.2e} photons/s/m²")

    # b) Scale the flux for a magnitude 20 star
    flux_ratio = 10**(-0.4 * star_mag)
    photon_rate_per_area_star = photon_rate_per_area_0 * flux_ratio

    print(f"\n  b) For a magnitude {star_mag} star:")
    print(f"    Photon flux rate (mag {star_mag}) = Rate(mag 0) * 10^(-0.4 * mag) = {photon_rate_per_area_0:.2e} * {flux_ratio:.1e} = {photon_rate_per_area_star:.2f} photons/s/m²")

    # c) Multiply by telescope area
    telescope_area_m2 = np.pi * (telescope_diameter_m / 2)**2
    total_photon_rate_r = photon_rate_per_area_star * telescope_area_m2
    
    print(f"\n  c) For a {telescope_diameter_m}m telescope:")
    print(f"    Telescope area = π * (D/2)² = π * ({telescope_diameter_m / 2})² = {telescope_area_m2:.3f} m²")
    print(f"    Total photon rate (R) = Flux * Area = {photon_rate_per_area_star:.2f} * {telescope_area_m2:.3f} = {total_photon_rate_r:.2f} photons/s")
    print("-" * 50)
    
    # 4. Calculate the exposure time
    exposure_time_s = total_photons_needed / total_photon_rate_r
    rounded_exposure_time = round(exposure_time_s)
    
    print("Step 4: Calculate Final Exposure Time (t)")
    print(f"Exposure Time (s) = Total Photons (N) / Total Photon Rate (R) = {total_photons_needed:.1f} / {total_photon_rate_r:.2f} = {exposure_time_s:.1f} s")
    print("\n" + "="*50)
    print(f"The required exposure time, rounded to the nearest integer, is: {rounded_exposure_time} seconds.")
    print("="*50)
    
    return rounded_exposure_time

if __name__ == '__main__':
    calculate_exposure_time()