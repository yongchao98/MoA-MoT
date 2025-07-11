import math

def calculate_transmission():
    """
    Calculates the transmission of characteristic X-rays through a Be window
    to determine the lightest detectable element in an EDX system.
    """
    # --- Constants ---
    # Beryllium (Be) window properties
    be_density_g_cm3 = 1.85  # Density (rho) in g/cm^3
    be_thickness_cm = 100 / 10000  # Thickness (x) in cm (100 um = 0.01 cm)

    # --- Element Data ---
    # Data for candidate elements: Name, K-alpha X-ray energy in keV, and
    # the mass attenuation coefficient (mu) of Be at that energy in cm^2/g.
    # Mu values are from NIST XCOM database for Be.
    elements = {
        'Na': {'energy_keV': 1.04, 'mu_be_cm2_g': 450.3},
        'Mg': {'energy_keV': 1.25, 'mu_be_cm2_g': 260.5},
        'Si': {'energy_keV': 1.74, 'mu_be_cm2_g': 110.1},
        'Ca': {'energy_keV': 3.69, 'mu_be_cm2_g': 16.48}
    }
    
    print("Calculating X-ray transmission through a 100 um Beryllium window...\n")
    
    # --- Calculation Loop ---
    for name, data in elements.items():
        mu = data['mu_be_cm2_g']
        
        # Beer-Lambert law exponent
        exponent = -mu * be_density_g_cm3 * be_thickness_cm
        
        # Transmission calculation
        transmission_fraction = math.exp(exponent)
        transmission_percent = transmission_fraction * 100
        
        print(f"--- For {name} (Energy = {data['energy_keV']} keV) ---")
        print(f"Equation: Transmission = exp(-mu * rho * x)")
        # Outputting each number in the final equation
        print(f"Calculation: exp(-{mu} cm^2/g * {be_density_g_cm3} g/cm^3 * {be_thickness_cm} cm)")
        print(f"Result: {transmission_percent:.2f}% transmission\n")
        
    print("--- Conclusion ---")
    print("An element is practically detectable if its X-ray transmission is at least a few percent.")
    print("Based on the calculations:")
    print("- Na and Mg have extremely low transmission (<1%) and are unlikely to be detected.")
    print("- Si has a transmission of ~13%, making it clearly detectable.")
    print("- Ca has a high transmission and is also easily detectable.")
    print("\nTherefore, Silicon (Si) is the lightest element among the choices that can be reliably seen.")

calculate_transmission()