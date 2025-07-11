import math

def estimate_photoelectron_density():
    """
    Estimates the density of photoelectrons in a cuvette of atomic hydrogen.
    """

    # 1. Constants and Input Parameters in SI units
    k_B = 1.380649e-23  # Boltzmann constant in J/K
    e = 1.602176634e-19   # Elementary charge in C
    
    # Given parameters
    T = 3000.0  # Temperature in Kelvin
    P_torr = 10.0  # Pressure in Torr
    I_w_cm2 = 10.0  # Intensity in W/cm^2

    # Convert inputs to SI units
    P_pa = P_torr * 133.322  # Pressure in Pascals (N/m^2)
    I_si = I_w_cm2 * 10000.0  # Intensity in W/m^2

    print("--- Input Parameters (SI Units) ---")
    print(f"Temperature (T): {T} K")
    print(f"Pressure (P): {P_pa:.2f} Pa")
    print(f"Intensity (I): {I_si:.1e} W/m^2")
    print("-" * 35)

    # 2. Photon Energy
    E_ion_eV = 13.6  # Ionization energy of Hydrogen in eV
    E_photon_eV = 2 * E_ion_eV  # Photon energy is 2 * E_ion = 27.2 eV
    E_photon_J = E_photon_eV * e  # Photon energy in Joules

    print("--- Photon Properties ---")
    print(f"Photon Energy (E_photon): {E_photon_J:.4e} J ({E_photon_eV} eV)")

    # 3. Photon Flux
    F = I_si / E_photon_J  # Photon flux in photons / (m^2 * s)
    print(f"Photon Flux (F): {F:.4e} photons/m^2/s")
    print("-" * 35)
    
    # 4. Hydrogen Atom Density (from Ideal Gas Law: P = n*k*T)
    n_H = P_pa / (k_B * T)
    print("--- Gas Properties ---")
    print(f"Initial Hydrogen Density (n_H): {n_H:.4e} atoms/m^3")
    
    # 5. Rate Coefficients
    # Photoionization cross-section for H at E=2*E_ion
    # sigma(E) ≈ sigma(E_ion) * (E_ion / E)^3
    # sigma(E_ion) ≈ 7.91e-18 cm^2 = 7.91e-22 m^2
    sigma_ion_at_threshold = 7.91e-22 # m^2
    sigma_ion = sigma_ion_at_threshold * (E_ion_eV / E_photon_eV)**3
    
    # Radiative recombination coefficient for H
    # alpha_rec(T) ≈ 2.7e-11 * T^(-1/2) cm^3/s for T in K
    # Convert to m^3/s: 2.7e-17 * T^(-1/2) m^3/s
    alpha_rec = 2.7e-17 * (T**(-0.5))

    print("--- Rate Coefficients ---")
    print(f"Photoionization Cross-section (σ_ion): {sigma_ion:.4e} m^2")
    print(f"Recombination Coefficient (α_rec): {alpha_rec:.4e} m^3/s")
    print("-" * 35)

    # 6. Solve the Rate Balance Equation
    # R_ionization = n_H * sigma_ion * F
    # R_recombination = n_e^2 * alpha_rec
    # Set R_ionization = R_recombination
    # n_e^2 = (n_H * sigma_ion * F) / alpha_rec
    
    numerator = n_H * sigma_ion * F
    n_e_squared = numerator / alpha_rec
    n_e = math.sqrt(n_e_squared)

    print("--- Calculation ---")
    print(f"The ionization rate per volume is {numerator:.4e} m^-3 s^-1.")
    print(f"The equation to solve is {n_H:.2e} * {sigma_ion:.2e} * {F:.2e} = {alpha_rec:.2e} * n_e^2.")
    print(f"This simplifies to n_e^2 = {numerator:.2e} / {alpha_rec:.2e} = {n_e_squared:.2e}.")
    print(f"The final calculated value for n_e is sqrt({n_e_squared:.2e}).")
    print("-" * 35)
    print("--- Final Answer ---")
    print(f"Estimated density of photoelectrons (n_e): {n_e:.4e} m^-3")

estimate_photoelectron_density()
<<<3.8499e+20>>>