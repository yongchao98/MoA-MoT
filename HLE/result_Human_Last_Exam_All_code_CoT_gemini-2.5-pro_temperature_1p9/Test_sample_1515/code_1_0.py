def print_nsvz_beta_function_sqcd(g_symbol, Nc, Nf):
    """
    Prints the NSVZ beta function formula for Supersymmetric QCD (SQCD).

    The NSVZ beta function is an exact result relating the beta function to the
    anomalous dimensions of the matter fields. Its derivation relies on properties
    like holomorphy, which must be preserved by the calculation scheme.

    The general formula is:
    β(g) = - [g³ / (16π²)] * [3T(adj) - Σ_i T(R_i)(1 - γ_i)] / [1 - (g²/(8π²))T(adj)]

    For SQCD with gauge group SU(Nc) and Nf flavors (Nf quarks and Nf anti-quarks):
    - T(adj) is the Dynkin index for the adjoint representation, which is Nc for SU(Nc).
    - The sum is over all 2*Nf matter chiral superfields. Each is in the fundamental
      representation with Dynkin index T(fund) = 1/2. So Σ T(R_i) = 2 * Nf * (1/2) = Nf.
    - We assume a common anomalous dimension γ for all matter fields.

    Args:
        g_symbol (str): The symbol to use for the gauge coupling (e.g., 'g').
        Nc (int): The number of colors.
        Nf (int): The number of flavors.
    """
    print("The exact NSVZ beta function for SU(Nc) SQCD with Nf flavors is given by:")
    print("β(g) = - [g³ / (16π²)] * [3*T(adj) - Nf * (1 - γ(g²))] / [1 - (g²/(8π²)) * T(adj)]")
    
    # For SU(Nc), the Dynkin index of the adjoint representation T(adj) is Nc.
    T_adj = Nc
    
    print(f"\nFor the specific case where the number of colors Nc = {Nc} and number of flavors Nf = {Nf}:")
    print(f"The required numbers for the equation are:")
    print(f"  - T(adj) = {T_adj}")
    print(f"  - Nf = {Nf}")
    
    # Constructing the string for the final equation with the numbers substituted.
    numerator_str = f"(3*{T_adj} - {Nf}*(1 - γ({g_symbol}²)))"
    denominator_str = f"(1 - {g_symbol}²/(8π²) * {T_adj})"
    
    print("\nThe final equation with these numbers is:")
    print(f"β({g_symbol}) = - ({g_symbol}³ / (16π²)) * {numerator_str} / {denominator_str}")


# --- Main Execution ---
# Define parameters for a sample supersymmetric gauge theory.
# The gauge coupling 'g' is treated as a variable.
gauge_coupling_symbol = 'g'
number_of_colors = 3      # Corresponds to gauge group SU(3)
number_of_flavors = 2     # Corresponds to 2 quark/anti-quark pairs

# Call the function to print the specific form of the NSVZ equation.
print_nsvz_beta_function_sqcd(gauge_coupling_symbol, number_of_colors, number_of_flavors)

<<<B>>>