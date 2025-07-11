import math

def calculate_transmission():
    """
    Calculates the X-ray transmission through a Be window for various elements
    and identifies the lightest one detectable.
    """
    # --- Setup Parameters ---
    # Density of Beryllium (g/cm^3)
    rho_be = 1.85
    # Thickness of the Beryllium window (cm)
    thickness_be_cm = 100e-4 # 100 micrometers = 0.01 cm

    print(f"Analyzing X-ray transmission through a {thickness_be_cm*1e4:.0f} µm Be window.\n")
    print(f"Using Beer-Lambert Law: T = exp(-μ * ρ * x)")
    print(f"Density of Be (ρ): {rho_be} g/cm^3")
    print(f"Window thickness (x): {thickness_be_cm} cm\n")
    print("-" * 50)

    # Data for elements:
    # - Name
    # - Z: Atomic Number
    # - Line: Characteristic X-ray line
    # - E (keV): Energy of the X-ray line
    # - mu_be (cm^2/g): Mass absorption coefficient of Be at energy E
    elements = [
        {'name': 'Na', 'Z': 11, 'line': 'Kα', 'E': 1.041, 'mu_be': 1234.0},
        {'name': 'Mg', 'Z': 12, 'line': 'Kα', 'E': 1.253, 'mu_be': 756.9},
        {'name': 'Si', 'Z': 14, 'line': 'Kα', 'E': 1.739, 'mu_be': 351.4},
        {'name': 'Ca', 'Z': 20, 'line': 'Kα', 'E': 3.691, 'mu_be': 58.33},
        {'name': 'W',  'Z': 74, 'line': 'Lα', 'E': 8.398, 'mu_be': 11.39} # Use L-alpha for heavy elements
    ]

    detectable_elements = []

    for el in sorted(elements, key=lambda i: i['Z']):
        mu = el['mu_be']
        exponent = mu * rho_be * thickness_be_cm
        transmission = math.exp(-exponent)
        transmission_percent = transmission * 100

        print(f"Element: {el['name']} (Z={el['Z']})")
        print(f"  X-ray Line: {el['line']} at {el['E']} keV")
        print(f"  μ of Be at this energy: {mu:.2f} cm^2/g")
        print(f"  Calculation: T = exp(-{mu:.2f} * {rho_be} * {thickness_be_cm})")
        print(f"  Transmission: {transmission_percent:.6f}%")

        # Set a detection threshold, e.g., 1%
        if transmission_percent > 1.0:
            print("  Result: Detectable\n")
            detectable_elements.append(el)
        else:
            print("  Result: Not practically detectable\n")
        print("-" * 50)

    if detectable_elements:
        # Find the element with the lowest Z among the detectable ones
        lightest_detectable = min(detectable_elements, key=lambda i: i['Z'])
        print(f"\nConclusion:")
        print(f"The lightest element from the list that can be seen is {lightest_detectable['name']}.")
    else:
        print("\nConclusion: None of the elements are detectable with this setup and a 1% threshold.")

if __name__ == "__main__":
    calculate_transmission()