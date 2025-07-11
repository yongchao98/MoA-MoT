import math

def calculate_nsvz_beta_pure_sym(g, N):
    """
    Calculates the NSVZ beta function for a pure N=1 Super Yang-Mills theory.

    Args:
        g (float): The gauge coupling constant.
        N (int): The number of colors for the SU(N) gauge group.
    """
    # The quadratic Casimir for the adjoint representation of SU(N) is T(G) = N.
    T_G = N

    # The numerator of the NSVZ formula for pure SYM is 3 * T(G).
    # The term with the anomalous dimension of matter is zero as there are no matter fields.
    numerator = 3 * T_G

    # The denominator involves a correction term.
    denominator = 1 - (g**2 * T_G) / (8 * math.pi**2)

    # The overall prefactor.
    prefactor = -g**3 / (16 * math.pi**2)

    # Calculate the beta function.
    beta_g = prefactor * (numerator / denominator)

    print("--- NSVZ Beta Function for Pure N=1 SU(N) SYM ---")
    print(f"For gauge group SU({N}) and coupling g = {g:.2f}:")
    print("\nCalculation steps:")
    print(f"1. Prefactor: -g^3 / (16*pi^2) = -({g:.2f}^3) / (16 * {math.pi:.4f}^2) = {-g**3:.4f} / {16*math.pi**2:.4f} = {prefactor:.4f}")
    print(f"2. Numerator term (3*T(G)): 3 * {T_G} = {numerator}")
    print(f"3. Denominator term (1 - g^2*T(G)/(8*pi^2)): 1 - ({g:.2f}^2 * {T_G}) / (8 * {math.pi:.4f}^2) = 1 - {g**2 * T_G:.4f} / {8*math.pi**2:.4f} = {denominator:.4f}")
    print(f"4. Final Beta function: beta(g) = {prefactor:.4f} * ({numerator} / {denominator:.4f})")
    print(f"\nResult:")
    print(f"beta(g) = {beta_g:.6f}")


# --- Main execution ---
# Example parameters for SU(3)
N_colors = 3
gauge_coupling = 1.2

# The crucial condition for this formula to hold is that the chosen regularization
# scheme must preserve the holomorphy of the theory. This allows the non-renormalization
# theorems to be fully leveraged.
print("The condition for the NSVZ beta function to match non-renormalization theorems is that the chosen quantum regularization scheme must preserve the holomorphy properties of the supersymmetric theory.")
print("-" * 60)

calculate_nsvz_beta_pure_sym(gauge_coupling, N_colors)
