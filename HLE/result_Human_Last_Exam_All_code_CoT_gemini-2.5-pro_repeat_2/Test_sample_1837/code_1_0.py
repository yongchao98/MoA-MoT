import numpy as np

def solve_star_measurement():
    """
    Analyzes the measurements of the star Pandora using Planck's Law.
    """
    # Physical constants
    h = 6.62607015e-34  # Planck's constant (J*s)
    c = 299792458      # Speed of light (m/s)
    k = 1.380649e-23   # Boltzmann constant (J/K)

    # Given measurements
    L_nm = 400.0       # Wavelength in nm
    L = L_nm * 1e-9    # Wavelength in meters
    T_meas = 9000.0    # Measured Temperature in K
    B_meas = 1.2e15    # Measured Spectral radiance in W/(m^2*sr*m)

    # Check for consistency using Planck's Law:
    # B(L, T) = (2 * h * c^2) / (L^5) * (1 / (exp(h*c / (L*k*T)) - 1))
    
    # Calculate theoretical spectral radiance with given T and L
    exponent_check = (h * c) / (L * k * T_meas)
    prefactor_check = (2 * h * c**2) / (L**5)
    B_calc = prefactor_check / (np.exp(exponent_check) - 1)

    # According to Wien's displacement law, the maximum possible spectral radiance
    # for a blackbody at T_meas=9000K occurs at a specific wavelength. Let's find
    # that maximum possible radiance.
    b_wien = 2.898e-3 # Wien's displacement constant in m*K
    L_max = b_wien / T_meas
    exponent_max = (h * c) / (L_max * k * T_meas)
    prefactor_max = (2 * h * c**2) / (L_max**5)
    B_max = prefactor_max / (np.exp(exponent_max) - 1)

    print(f"The measured spectral radiance is B = {B_meas:.2e} W/(m^2*sr*m).")
    print(f"The calculated spectral radiance for T=9000 K and L=400 nm is B_calc = {B_calc:.2e} W/(m^2*sr*m).")
    print(f"The maximum possible radiance for a 9000 K blackbody is {B_max:.2e} W/(m^2*sr*m).")
    print("The measured radiance is physically impossible for the given temperature.")
    print("This suggests the temperature measurement is incorrect, as 9000 K is also unusually low for a DB-class white dwarf.")
    
    # Assume L and B are correct and calculate the expected temperature T.
    # T = (h*c) / (L*k * ln( (2*h*c^2 / (L^5 * B_meas)) + 1 ))
    print("\nCalculating the expected temperature assuming L and B are correct...")
    print("Equation: T = (h * c) / (L * k * ln( (2 * h * c^2) / (L^5 * B) + 1 ))")
    print("\nUsing the values:")
    print(f"h = {h:.4e} J*s")
    print(f"c = {c:.4e} m/s")
    print(f"k = {k:.4e} J/K")
    print(f"L = {L:.4e} m")
    print(f"B = {B_meas:.2e} W/(m^2*sr*m)")
    
    term_in_ln = (prefactor_check / B_meas) + 1
    T_calc = (h * c) / (L * k * np.log(term_in_ln))
    
    # Final answer
    X = 'T'
    v = int(round(T_calc))
    
    print(f"\nThe corrected temperature is {T_calc:.0f} K.")
    print("This value is within the typical range for a DB-class white dwarf (12,000 K - 30,000 K).")
    print("Therefore, the quantity most likely having a measurement error is T.")
    
    print(f"\n<<<{X}{v}>>>")

solve_star_measurement()