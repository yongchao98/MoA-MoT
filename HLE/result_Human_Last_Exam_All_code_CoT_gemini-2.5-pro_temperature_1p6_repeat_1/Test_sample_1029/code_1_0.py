def compute_poynting_vector_expression():
    """
    This function prints the symbolic expression for the Poynting vector
    for a long, charged, insulating cylindrical rod moving along its axis
    within an external uniform electric field.

    The derivation considers two regions:
    1. Inside the rod (where radial distance r < Radius R)
    2. Outside the rod (where radial distance r >= Radius R)

    The final Poynting vector S is presented with its radial (r_hat) and
    axial (k_hat) components.

    Variables in the equations:
    E:         Magnitude of the external electric field.
    rho:       Uniform volume charge density of the rod.
    v:         Speed of the rod along its axis.
    R:         Radius of the rod.
    r:         Radial distance from the rod's axis.
    epsilon_0: Permittivity of free space.
    """

    print("The Poynting vector S = S_r * r_hat + S_k * k_hat is computed for two regions.")
    print("="*70)

    # --- Region 1: Inside the rod (r < R) ---
    print("\nRegion 1: Inside the rod (for r < R)\n")
    print("-"*70)

    # Radial component S_r
    print("The radial component (S_r) of the Poynting vector is:")
    # The formula is S_r = (E * rho * v * r) / 2
    print("S_r = (1/2) * E * rho * v * r")
    print()

    # Axial component S_k
    print("The axial component (S_k) of the Poynting vector is:")
    # The formula is S_k = (rho^2 * v * r^2) / (4 * epsilon_0)
    print("S_k = (1/4) * (rho**2 * v * r**2) / epsilon_0")
    print()

    # Full vector expression
    print("The complete Poynting vector inside the rod is:")
    print("S(r) = [(1/2) * E * rho * v * r] r_hat  +  [(1/4) * (rho**2 * v * r**2) / epsilon_0] k_hat")
    print("="*70)

    # --- Region 2: Outside the rod (r >= R) ---
    print("\nRegion 2: Outside the rod (for r >= R)\n")
    print("-"*70)

    # Radial component S_r
    print("The radial component (S_r) of the Poynting vector is:")
    # The formula is S_r = (E * rho * v * R^2) / (2 * r)
    print("S_r = (1/2) * (E * rho * v * R**2) / r")
    print()

    # Axial component S_k
    print("The axial component (S_k) of the Poynting vector is:")
    # The formula is S_k = (rho^2 * v * R^4) / (4 * epsilon_0 * r^2)
    print("S_k = (1/4) * (rho**2 * v * R**4) / (epsilon_0 * r**2)")
    print()

    # Full vector expression
    print("The complete Poynting vector outside the rod is:")
    print("S(r) = [(1/2) * (E * rho * v * R**2) / r] r_hat  +  [(1/4) * (rho**2 * v * R**4) / (epsilon_0 * r**2)] k_hat")
    print("="*70)


if __name__ == "__main__":
    compute_poynting_vector_expression()