import numpy as np

def solve_bose_equilibrium():
    """
    Calculates and prints the equilibrium equations for entropy and mean energy
    for the Bose case of light quanta, justified by large deviation theorems.
    """

    # --- Introduction and Theoretical Justification ---
    print("--- Equilibrium Relations for Light Quanta (Bose-Einstein Statistics) ---")
    print("""
The equilibrium state of a physical system is its most probable macroscopic state.
Large deviation theorems, like Sanov's theorem, provide a formal justification
for the Principle of Maximum Entropy, which states that the equilibrium distribution
is the one that maximizes entropy (S) subject to physical constraints (e.g., mean energy U).

For bosons like light quanta, this procedure yields the following key relationships:
    """)
    print("Let:")
    print("  U       = Equilibrium Mean Energy of a single mode")
    print("  S       = Equilibrium Entropy of a single mode")
    print("  epsilon = Energy of a single quantum")
    print("  k_B     = Boltzmann constant")
    print("  T       = Temperature")
    print("  beta    = 1 / (k_B * T) is the inverse temperature\n")

    # --- Present the Derived Equations ---
    u_equation_str = "U = epsilon / (exp(beta * epsilon) - 1)"
    print("1. Equilibrium Mean Energy (U), known as the Planck distribution:")
    print(f"   Equation: {u_equation_str}\n")

    s_equation_str = "S / k_B = (1 + U/epsilon) * log(1 + U/epsilon) - (U/epsilon) * log(U/epsilon)"
    print("2. Equilibrium Entropy (S) as a function of Mean Energy (U):")
    print(f"   Equation: {s_equation_str}\n")


    # --- Provide a Numerical Example ---
    print("--- Numerical Example Calculation ---")
    
    # Set dimensionless values for the constants for the example.
    k_B = 1.0      # Boltzmann constant
    epsilon = 1.0  # Energy of a quantum
    T = 1.0        # Temperature
    beta = 1.0 / (k_B * T)

    print(f"We assume a system where epsilon = {epsilon:.1f}, k_B = {k_B:.1f}, and T = {T:.1f} (so beta = {beta:.1f}).\n")

    # Calculate U with step-by-step printing
    print("--> Calculating Mean Energy (U):")
    U = epsilon / (np.exp(beta * epsilon) - 1)
    
    # Show each number in the equation
    print(f"    U = {epsilon:.4f} / (exp({beta:.4f} * {epsilon:.4f}) - 1)")
    exp_val = np.exp(beta * epsilon)
    print(f"    U = {epsilon:.4f} / ({exp_val:.4f} - 1)")
    denominator_u = exp_val - 1
    print(f"    U = {epsilon:.4f} / {denominator_u:.4f}")
    print(f"    U = {U:.4f}\n")

    # Calculate S with step-by-step printing
    print("--> Calculating Entropy (S) from the calculated U:")
    
    n_mean = U / epsilon # Mean number of quanta
    
    # Handle the n_mean=0 case to avoid log(0)
    if n_mean < 1e-9:
        S = 0.0
        print("    U is near zero, so S is zero.")
    else:
        # Show each number in the equation
        print(f"    S / {k_B:.4f} = (1 + {n_mean:.4f}) * log(1 + {n_mean:.4f}) - ({n_mean:.4f}) * log({n_mean:.4f})")
        
        # Calculate intermediate terms
        log_term1_arg = 1 + n_mean
        log_term1_val = np.log(log_term1_arg)
        log_term2_arg = n_mean
        log_term2_val = np.log(log_term2_arg)

        print(f"    S / {k_B:.4f} = ({log_term1_arg:.4f}) * log({log_term1_arg:.4f}) - ({log_term2_arg:.4f}) * log({log_term2_arg:.4f})")
        print(f"    S / {k_B:.4f} = ({log_term1_arg:.4f}) * {log_term1_val:.4f} - ({log_term2_arg:.4f}) * ({log_term2_val:.4f})")
        
        term1 = log_term1_arg * log_term1_val
        term2 = log_term2_arg * log_term2_val
        print(f"    S / {k_B:.4f} = {term1:.4f} - ({term2:.4f})")
        
        S_over_kB = term1 - term2
        print(f"    S / {k_B:.4f} = {S_over_kB:.4f}")
        
        S = k_B * S_over_kB
        print(f"    S = {S:.4f}")
    
    print("\n--- Final Equilibrium Values for the Example ---")
    print(f"For a system with beta = {beta:.1f}, the equilibrium mean energy is {U:.4f} and entropy is {S:.4f}.")
    print(f"<<<U={U:.4f}, S={S:.4f}>>>")


if __name__ == '__main__':
    solve_bose_equilibrium()