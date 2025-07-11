def solve_partition_function():
    """
    This function prints a step-by-step derivation of the partition function Z
    for a system with Hamiltonian H = -μN in the grand canonical ensemble.
    """
    
    print("Derivation of the Partition Function Z\n")
    
    # Step 1: Definition of the Grand Canonical Partition Function
    print("Step 1: The grand canonical partition function Z is defined as a trace over the system's states:")
    print("Z = Tr[exp(-β * (H - μ*N))]")
    print("Here, β is the inverse temperature, μ is the chemical potential, H is the Hamiltonian, and N is the number operator.\n")
    
    # Step 2: Substitute the given Hamiltonian
    print("Step 2: Substitute the given Hamiltonian H = -μ*N into the formula.")
    print("Z = Tr[exp(-β * ((-μ*N) - μ*N))]")
    print("This simplifies the term in the exponent:\n")
    print("-β * (-μ*N - μ*N) = -β * (-2*μ*N) = 2*β*μ*N\n")
    
    # Step 3: Simplified expression for Z
    print("Step 3: The expression for Z becomes:")
    print("Z = Tr[exp(2*β*μ*N)]\n")

    # Step 4: Evaluate the Trace
    print("Step 4: The trace is evaluated by summing over a basis of states. We use the number basis |n>, where n is the number of particles.")
    print("N|n> = n|n>")
    print("So, Tr[exp(2*β*μ*N)] = Σ_n exp(2*β*μ*n)")
    print("The allowed values for the occupation number 'n' depend on the particle statistics.\n")
    
    print("--------------------------------------------------\n")

    # Case A: Fermions
    print("Case A: Fermionic System (for a single quantum state)\n")
    print("For fermions, the Pauli exclusion principle applies, so a single state can be either empty (n=0) or occupied by one particle (n=1).\n")
    print("Z_F = Σ_{n=0,1} exp(2 * β * μ * n)")
    print("   = exp(2 * β * μ * 0) + exp(2 * β * μ * 1)")
    print("   = 1 + exp(2 * β * μ)\n")
    
    # Case B: Bosons
    print("Case B: Bosonic System (for a single quantum state)\n")
    print("For bosons, any number of particles can occupy a single state (n = 0, 1, 2, ...).\n")
    print("Z_B = Σ_{n=0 to ∞} exp(2 * β * μ * n)")
    print("This is a geometric series: Σ_{n=0 to ∞} [exp(2*β*μ)]^n")
    print("The series converges to 1 / (1 - x) if |x| < 1. This requires exp(2*β*μ) < 1, which means μ must be negative.\n")
    print("Z_B = 1 / (1 - exp(2 * β * μ))  (for μ < 0)\n")

    print("--------------------------------------------------\n")

    # Path Integral Connection
    print("Connection to Path Integrals:\n")
    print("The path integral method calculates the same partition function Z by evaluating a functional integral.")
    print("For this system, the path integral action S is:")
    print("S = ∫dτ d^d x [ ψ*(∂_τ - 2μ)ψ ]")
    print("Calculating Z = ∫D[ψ*,ψ] exp(-S) for fermions and bosons yields the exact same results Z_F and Z_B found above.\n")

    print("Final Answer Summary (for a single mode):")
    print("  - For Fermions: Z = 1 + exp(2 * β * μ)")
    print("  - For Bosons:   Z = 1 / (1 - exp(2 * β * μ))")


if __name__ == "__main__":
    solve_partition_function()
