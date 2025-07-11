import numpy as np

def calculate_nsvz_beta_pure_sym(g, Nc):
    """
    Calculates the NSVZ beta function for a pure N=1 Super-Yang-Mills theory.

    The NSVZ formula for a pure N=1 SYM theory is:
    β(g) = - [g^3 / (16 * π^2)] * [3 * T(adj)] / [1 - (T(adj) * g^2) / (8 * π^2)]

    Args:
        g (float): The gauge coupling constant.
        Nc (int): The number of colors for the SU(Nc) gauge group.

    Returns:
        float: The value of the beta function.
    """
    # For SU(Nc), the Dynkin index of the adjoint representation is T(adj) = Nc.
    T_adj = float(Nc)

    # Calculate the pre-factor
    prefactor = -g**3 / (16 * np.pi**2)

    # Calculate the numerator term (no matter fields, so T(R)=0)
    numerator = 3 * T_adj

    # Calculate the denominator term
    denominator = 1 - (T_adj * g**2) / (8 * np.pi**2)

    # Check for potential division by zero (Landau pole)
    if abs(denominator) < 1e-9:
        print("Warning: Denominator is close to zero; approaching a Landau pole.")
        return float('inf')

    # Calculate the beta function
    beta_g = prefactor * (numerator / denominator)

    # --- Output the equation with the calculated numbers ---
    print(f"Calculating the NSVZ beta function for SU({Nc}) pure N=1 SYM with g = {g:.2f}")
    print(f"The Dynkin index for the adjoint representation, T(adj), is: {T_adj:.2f}")
    print("\nThe NSVZ Beta Function Formula:")
    print("β(g) = -[g³ / (16π²)] * [3 * T(adj)] / [1 - T(adj)g² / (8π²)]")
    print("\nPlugging in the values:")
    
    # We break down the calculation for clarity.
    num_calc = 3 * T_adj
    den_calc_term = T_adj * g**2 / (8 * np.pi**2)
    den_calc = 1 - den_calc_term
    pre_calc = -g**3 / (16 * np.pi**2)

    # Final formatted output of the equation
    print(f"β({g:.2f}) = -[{g**3:.2f} / (16 * {np.pi**2:.2f})] * [{num_calc:.2f}] / [1 - ({T_adj:.2f} * {g**2:.2f}) / (8 * {np.pi**2:.2f})]")
    print(f"β({g:.2f}) = [{pre_calc:.4f}] * [{num_calc:.2f}] / [{den_calc:.4f}]")
    print(f"β({g:.2f}) = {beta_g:.6f}")
    
    return beta_g

# --- Main execution ---
if __name__ == "__main__":
    # Define parameters for the calculation
    gauge_coupling = 1.5  # A value for the coupling constant
    number_of_colors = 3  # For gauge group SU(3), like in QCD
    
    calculate_nsvz_beta_pure_sym(g=gauge_coupling, Nc=number_of_colors)