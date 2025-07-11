def print_fermionic_partition_function():
    """
    Prints the formula for the fermionic partition function Z in the
    imaginary time path integral formalism.
    """
    
    # Using Unicode characters for mathematical symbols to make the formula clear.
    integral_symbol = "\u222B"
    d_symbol = "D"
    psi_bar_symbol = "\u03C8\u0304"  # psi with a bar over it
    psi_symbol = "\u03C8"
    exp_symbol = "exp"
    beta_symbol = "\u03B2"
    tau_symbol = "\u03C4"
    partial_symbol = "\u2202"
    mu_symbol = "\u03BC"

    print("The formula for the fermionic partition function Z is given by a functional path integral over Grassmann fields.")
    print("-" * 80)
    
    # Step 1: Print the overall structure of the partition function Z
    print("1. The partition function Z is defined as the integral over all field configurations:")
    z_definition = f"Z = {integral_symbol} {d_symbol}{psi_bar_symbol} {d_symbol}{psi_symbol} {exp_symbol}[-S_E]"
    print(f"\n  {z_definition}\n")
    print(f"   Where:")
    print(f"   - '{integral_symbol} {d_symbol}{psi_bar_symbol} {d_symbol}{psi_symbol}' denotes the functional integral over all anti-commuting Grassmann fields {psi_symbol} and {psi_bar_symbol}.")
    print(f"   - 'S_E' is the Euclidean action.\n")

    # Step 2: Define the Euclidean Action S_E
    print("2. The Euclidean Action, S_E, for a general non-interacting system is:")
    s_e_definition = f"S_E = {integral_symbol}\u2080\u02E3 d{tau_symbol} {integral_symbol} d\u00B3x [ {psi_bar_symbol}({tau_symbol}, x) ({partial_symbol}_{tau_symbol} + H - {mu_symbol}) {psi_symbol}({tau_symbol}, x) ]"
    print(f"\n  {s_e_definition}\n")
    print(f"   Where:")
    print(f"   - '{integral_symbol}\u2080\u02E3' is the integral over imaginary time '{tau_symbol}' from 0 to '{beta_symbol}'.")
    print(f"   - '{beta_symbol} = 1 / (k_B T)' is the inverse temperature.")
    print(f"   - '{integral_symbol} d\u00B3x' is the integral over spatial dimensions.")
    print(f"   - '{psi_bar_symbol}({tau_symbol}, x)' and '{psi_symbol}({tau_symbol}, x)' are the Grassmann fields.")
    print(f"   - '{partial_symbol}_{tau_symbol}' is the partial derivative with respect to imaginary time.")
    print(f"   - 'H' is the single-particle Hamiltonian operator (e.g., -∇²/2m).")
    print(f"   - '{mu_symbol}' is the chemical potential.\n")

    # Step 3: Specify the boundary conditions
    print("3. The path integral for fermions requires anti-periodic boundary conditions in the imaginary time direction:")
    boundary_condition = f"{psi_symbol}({tau_symbol} = {beta_symbol}, x) = -{psi_symbol}({tau_symbol} = 0, x)"
    print(f"\n  {boundary_condition}\n")
    print("   This condition is a direct consequence of the fermionic anti-commutation relations when evaluating the trace in the definition Z = Tr(exp(-βH)).")
    print("-" * 80)


if __name__ == "__main__":
    print_fermionic_partition_function()