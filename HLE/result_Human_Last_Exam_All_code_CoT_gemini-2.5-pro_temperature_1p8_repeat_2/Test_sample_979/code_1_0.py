def display_magnetic_field_solution():
    """
    Prints the symbolic solution for the magnetic field H(r, theta)
    for a spherical shell with surface current.
    """
    print("Based on the derivation, the magnetic field H(r, theta) is:")
    print("-" * 60)

    # --- Inside the sphere ---
    print("For the region inside the sphere (0 < r < R):")
    print("H(r, theta) = ( (2 * mu_0 / mu) / (1 + 2 * mu_0 / mu) ) * K_0 * z_hat")
    print("\nThis describes a uniform magnetic field in the z-direction.")
    
    print("-" * 60)

    # --- Outside the sphere ---
    print("For the region outside the sphere (r > R):")
    print("H(r, theta) = ( K_0 / (1 + 2 * mu_0 / mu) ) * (R^3 / r^3) * (2*cos(theta)*r_hat + sin(theta)*theta_hat)")
    print("\nThis describes a magnetic dipole field.")

    print("-" * 60)
    print("This corresponds to Choice E.")
    print("Where:")
    print("  K_0: Current amplitude")
    print("  mu: Magnetic permeability inside the sphere")
    print("  mu_0: Magnetic permeability of free space")
    print("  R: Radius of the sphere")
    print("  (r, theta): Spherical coordinates")
    print("  z_hat, r_hat, theta_hat: Unit vectors")

if __name__ == '__main__':
    display_magnetic_field_solution()