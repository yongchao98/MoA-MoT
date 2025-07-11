import math

def solve_edx_problem():
    """
    Calculates the X-ray transmission for various elements through a Be window
    to determine the lightest detectable element in an EDX system.
    """
    # --- System Parameters ---
    # Density of Beryllium (g/cm^3)
    rho_be = 1.85
    # Thickness of the Be window (cm)
    x_be = 100e-4  # 100 micrometers = 0.01 cm

    # --- Element Data ---
    # Data includes:
    # 'Z': Atomic number
    # 'E_ka': K-alpha X-ray energy (keV)
    # 'mu_be_at_E_ka': Mass attenuation coefficient of Be at the element's K-alpha energy (cm^2/g)
    elements = {
        'Na': {'Z': 11, 'E_ka': 1.041, 'mu_be_at_E_ka': 680},
        'Mg': {'Z': 12, 'E_ka': 1.254, 'mu_be_at_E_ka': 395},
        'Si': {'Z': 14, 'E_ka': 1.740, 'mu_be_at_E_ka': 160},
        'Ca': {'Z': 20, 'E_ka': 3.692, 'mu_be_at_E_ka': 21.5},
    }

    print("Calculating X-ray transmission through a 100 µm Be window.")
    print("Formula: Transmission = exp(-μ * ρ * x)\n")
    print(f"Parameters: ρ(Be) = {rho_be} g/cm³, x(window) = {x_be*1e4} µm = {x_be} cm\n")

    results = {}
    # Sort elements by atomic number to evaluate from lightest to heaviest
    sorted_elements = sorted(elements.items(), key=lambda item: item[1]['Z'])

    for name, data in sorted_elements:
        mu = data['mu_be_at_E_ka']
        exponent = -mu * rho_be * x_be
        transmission = math.exp(exponent)
        results[name] = transmission

        print(f"--- Element: {name} (Z={data['Z']}) ---")
        print(f"Kα Energy = {data['E_ka']} keV")
        print(f"μ in Be at this energy = {mu} cm²/g")
        print(f"Calculation: exp(-{mu} * {rho_be} * {x_be})")
        print(f"Resulting Transmission = {transmission:.6f} or {transmission:.4%}")
        print("-" * 25 + "\n")

    # Determine the lightest element with transmission > 1%
    lightest_detectable = None
    for name, data in sorted_elements:
        if results[name] > 0.01:  # 1% transmission threshold
            lightest_detectable = name
            break
            
    print("Conclusion:")
    print("A common threshold for detectability is a transmission greater than 1%.")
    print(f"Based on the calculations, Na and Mg have transmissions far below 1%.")
    print(f"Si has a transmission of {results['Si']:.4%}, which is clearly detectable.")
    print(f"Therefore, Si is the lightest element among the choices that can be reliably seen.")


solve_edx_problem()
<<<E>>>