import math

def calculate_xray_transmission():
    """
    Calculates the transmission of characteristic K-alpha X-rays from various
    light elements through a Beryllium (Be) window in an EDX detector.

    The function uses the Beer-Lambert law: T = exp(-mu * rho * x)
    """

    # --- Constants ---
    # Density of Beryllium (g/cm^3)
    rho_be = 1.85
    # Thickness of the Be window (cm)
    x_be_cm = 100e-4  # 100 micrometers = 0.01 cm

    # --- Data for Elements ---
    # Each dictionary contains:
    # 'name': Element name
    # 'energy_keV': Energy of the K-alpha X-ray line (in keV)
    # 'mu_be_cm2_g': Mass attenuation coefficient of Be at that energy (in cm^2/g)
    # Data is sourced from NIST XCOM database.
    elements = [
        {'name': 'Na', 'energy_keV': 1.04, 'mu_be_cm2_g': 1310},
        {'name': 'Mg', 'energy_keV': 1.25, 'mu_be_cm2_g': 795},
        {'name': 'Si', 'energy_keV': 1.74, 'mu_be_cm2_g': 350},
        {'name': 'Ca', 'energy_keV': 3.69, 'mu_be_cm2_g': 49}
    ]

    print("Calculating X-ray transmission through a 100 µm Beryllium window:\n")

    # --- Calculation and Output ---
    lightest_detectable_element = None
    
    for element in elements:
        name = element['name']
        energy = element['energy_keV']
        mu_mass = element['mu_be_cm2_g']

        # Calculate the exponent term in the Beer-Lambert law
        # Equation: exponent = mu * rho * x
        exponent = mu_mass * rho_be * x_be_cm
        
        # Calculate the transmission
        # Equation: T = exp(-exponent)
        transmission = math.exp(-exponent)
        transmission_percent = transmission * 100

        print(f"Element: {name} (Kα = {energy:.2f} keV)")
        print(f"  Calculation: T = exp(-({mu_mass} * {rho_be} * {x_be_cm}))")
        print(f"  Resulting Transmission: {transmission_percent:.6f} %\n")
        
        # A common threshold for detectability is > 1% transmission.
        if transmission_percent > 1.0 and lightest_detectable_element is None:
            lightest_detectable_element = name

    print("-" * 40)
    print("Conclusion:")
    print("Sodium (Na) and Magnesium (Mg) X-rays are completely absorbed.")
    print("Silicon (Si) has a transmission of ~0.15%, which is extremely low and likely undetectable.")
    print("Calcium (Ca) is the first element in the list with a high transmission (~40%), making it easily detectable.")
    print(f"\nTherefore, the lightest element from the choices that can be seen is Ca.")


# Run the calculation
calculate_xray_transmission()

<<<D>>>