import sys

def solve_partition_function():
    """
    This function prints the path integral expression for the grand canonical partition function.
    """
    # Using Unicode characters for mathematical symbols to improve readability.
    integral_symbol = "\u222B"
    beta_symbol = "\u03B2"
    tau_symbol = "\u03C4"
    mu_symbol = "\u03BC"
    psi_star_symbol = "\u03C8*"
    psi_symbol = "\u03C8"
    partial_tau_symbol = "\u2202_" + tau_symbol
    hbar_symbol = "\u0127"
    nabla_squared_symbol = "\u2207\u00B2"
    
    # Define the components of the equations.
    # We describe the fields as functions of space (x) and imaginary time (tau).
    field_vars = f"({', '.join(['x', tau_symbol])})"
    psi_star_full = f"{psi_star_symbol}{field_vars}"
    psi_full = f"{psi_symbol}{field_vars}"
    
    # The partition function Z is a functional integral over the fields.
    path_integral_measure = f"{integral_symbol} D[{psi_star_symbol}] D[{psi_symbol}]"
    partition_function_z = f"Z = {path_integral_measure} exp(-S_E)"
    
    # The Euclidean action S_E is an integral of the Euclidean Lagrangian L_E.
    # For a general single-particle Hamiltonian operator H_op, the action is:
    integrals_s_e = f"{integral_symbol}\u2080^{beta_symbol} d{tau_symbol} {integral_symbol} d\u00B3x"
    lagrangian_content = f"{psi_star_full} [ {partial_tau_symbol} + H_op - {mu_symbol} ] {psi_full}"
    euclidean_action_s_e = f"S_E = {integrals_s_e} ( {lagrangian_content} )"

    # --- Printing the final result ---
    print("The path integral representation for the grand canonical partition function Z is given by the following equations:")
    print("-" * 80)

    print("\n1. Partition Function (Z):")
    print(f"   {partition_function_z}")
    
    print("\n2. Euclidean Action (S_E):")
    print(f"   {euclidean_action_s_e}")
    
    print("-" * 80)
    print("\nExplanation of the terms:")
    print(f"  Z: The grand canonical partition function.")
    print(f"  {path_integral_measure}: The path integral over all possible configurations of the fields {psi_symbol} and {psi_star_symbol}.")
    print(f"  S_E: The Euclidean action of the system.")
    print(f"  {beta_symbol} (beta): Inverse temperature, equal to 1/(k_B T).")
    print(f"  {tau_symbol} (tau): Imaginary time, which is integrated over the interval [0, {beta_symbol}].")
    print(f"  x: The spatial coordinates, integrated over the volume of the system.")
    print(f"  {psi_symbol}, {psi_star_symbol}: The matter fields representing the particles (e.g., complex scalar fields for bosons).")
    print(f"  {partial_tau_symbol}: The partial derivative with respect to imaginary time.")
    print(f"  H_op: The single-particle Hamiltonian operator. For a non-relativistic free particle, H_op = -({hbar_symbol}\u00B2/2m){nabla_squared_symbol}.")
    print(f"  {mu_symbol} (mu): The chemical potential of the system.")
    print("-" * 80)

# Execute the function to display the result.
solve_partition_function()