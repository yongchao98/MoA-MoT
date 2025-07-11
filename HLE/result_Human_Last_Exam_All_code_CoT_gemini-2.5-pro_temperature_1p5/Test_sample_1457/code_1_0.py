import math

def calculate_transmission():
    """
    Calculates the transmission of X-rays through a Beryllium window
    to determine the lightest detectable element in an EDX system.
    """
    # --- System Parameters ---
    # Density of Beryllium (g/cm^3)
    rho_be = 1.85
    # Thickness of the Be window (cm)
    t_be = 0.01 # 100 micrometers = 0.01 cm

    # --- Element Data ---
    # K-alpha X-ray energies (keV) and mass attenuation coefficients for Be at those energies (cm^2/g)
    # Source for coefficients: NIST XCOM database
    elements = {
        'Na': {'Z': 11, 'E_ka': 1.041, 'mu_rho_be': 615.1},
        'Mg': {'Z': 12, 'E_ka': 1.253, 'mu_rho_be': 358.3},
        'Si': {'Z': 14, 'E_ka': 1.740, 'mu_rho_be': 151.0},
        'Ca': {'Z': 20, 'E_ka': 3.691, 'mu_rho_be': 26.39}
    }

    print("Calculating X-ray transmission through a 100 um Beryllium window...\n")
    print(f"Using formula: Transmission = exp(-(mu/rho) * rho * t)")
    print(f"With Beryllium density (rho) = {rho_be} g/cm^3 and window thickness (t) = {t_be} cm.\n")
    
    lightest_detectable = None
    min_detectable_z = float('inf')

    # Sort elements by atomic number (Z) to find the lightest
    sorted_elements = sorted(elements.items(), key=lambda item: item[1]['Z'])

    for name, data in sorted_elements:
        mu_rho = data['mu_rho_be']
        exponent_value = -mu_rho * rho_be * t_be
        transmission = math.exp(exponent_value)
        transmission_percent = transmission * 100

        print(f"--- Element: {name} (Z={data['Z']}) ---")
        print(f"Kα Energy = {data['E_ka']} keV")
        print(f"Be Mass Attenuation Coefficient (mu/rho) = {mu_rho} cm^2/g")
        
        # Display the full calculation as requested
        print(f"Transmission = exp(-({mu_rho} * {rho_be} * {t_be}))")
        print(f"Transmission = exp({exponent_value:.3f})")
        print(f"Result = {transmission:.6f} or {transmission_percent:.4f}%\n")

        # A practical threshold for detection is typically >1% transmission.
        if transmission_percent > 1.0 and data['Z'] < min_detectable_z:
            min_detectable_z = data['Z']
            lightest_detectable = name
    
    print("--- Conclusion ---")
    print("Sodium (Na) and Magnesium (Mg) have transmission percentages well below 1%, making them practically undetectable with this thick Be window.")
    print(f"Silicon (Si) is the first element in the list (ordered by atomic number) with a significant transmission ({elements['Si']['transmission']:.2f}%).".replace("'",'"'))
    print(f"Therefore, Silicon is the lightest element among the choices that can be reliably seen.")


# Create a mock 'transmission' key for the final print statement
elements_for_print = {'Si': {'transmission': 100 * math.exp(-151.0 * 1.85 * 0.01)}}
calculate_transmission.__globals__['elements']['Si']['transmission'] = elements_for_print['Si']['transmission']

if __name__ == '__main__':
    calculate_transmission()
