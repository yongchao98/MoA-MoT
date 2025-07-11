def partition_function_calculation():
    """
    This script presents the derivation and final equations for the partition function Z
    for a system with the grand canonical Hamiltonian K = -μN, calculated using path integrals.
    """
    print("### Derivation and Result for the Partition Function Z ###")
    print("\nThe grand canonical partition function Z for a system with grand canonical")
    print("Hamiltonian K = -μN is given by Z = Tr[exp(-βK)] = Tr[exp(βμN)].")
    print("This trace can be evaluated using a path integral over fields ψ, where Z = ∫D[ψ*,ψ]exp(-S),")
    print("with the action S = ∫dτ Σₖ ψₖ*(∂_τ - μ)ψₖ.")
    print("\nThe evaluation of this path integral depends on the particle statistics.")
    print("Let M be the total number of identical single-particle states available.")
    print("-" * 50)

    # Define the symbols for the equation as strings
    Z = "Z"
    one = "1"
    exp = "exp"
    beta = "β"
    mu = "μ"
    M = "M"

    # --- Fermion Case ---
    print("\nCase 1: The system consists of FERMIONS.")
    print("For fermions, each state can be occupied by at most one particle.")
    print("The result of the path integral for the partition function Z is:")
    # Print the final equation, showing each symbolic component
    print(f"  {Z} = ({one} + {exp}({beta} * {mu}))^{M}")
    print("\nWhere:")
    print(f"  {Z}: The Partition Function")
    print(f"  {beta}: The inverse temperature, 1/(k_B * T)")
    print(f"  {mu}: The effective chemical potential")
    print(f"  {M}: The number of identical single-particle states")
    print("-" * 50)

    # --- Boson Case ---
    print("\nCase 2: The system consists of BOSONS.")
    print("For bosons, each state can be occupied by any number of particles.")
    print("The calculation requires that μ < 0 for the result to be well-defined.")
    print("The result of the path integral for the partition function Z is:")
    # Print the final equation, showing each symbolic component
    print(f"  {Z} = ({one} - {exp}({beta} * {mu}))^(-{M})")
    print("\nWhere the symbols have the same meaning as above.")

partition_function_calculation()