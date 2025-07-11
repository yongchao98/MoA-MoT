def display_fermionic_partition_function():
    """
    Prints the formula for the fermionic partition function Z
    in the imaginary time path integral representation.
    """
    print("The formula for the fermionic partition function (Z) in the imaginary time path integral formalism is as follows:")
    print("=" * 80)

    # 1. The main partition function formula
    print("\n1. The Partition Function (Z):")
    print("   The partition function is a functional integral over all anti-commuting Grassmann fields ψ and ψ̄.")
    print("\n   Z = ∫ D[ψ̄]D[ψ] * e^(-S[ψ̄, ψ])\n")

    # 2. The Euclidean Action (S)
    print("2. The Euclidean Action (S):")
    print("   The weight of each path is determined by the Euclidean action, which is the integral of the Euclidean Lagrangian (L_E).")
    print("\n   S[ψ̄, ψ] = ∫[from τ=0 to β] dτ * ∫ dᵈx * L_E(ψ̄, ψ)\n")

    # 3. The Boundary Conditions
    print("3. Fermionic Boundary Condition:")
    print("   A key feature for fermions is the anti-periodic boundary condition in imaginary time (τ), which arises from the trace operation.")
    print("   The fields at the beginning and end of the imaginary time interval are related by:")
    print("\n   ψ(τ=β) = -1 * ψ(τ=0)\n")
    print("   Note the numbers β, -1, and 0 which define this crucial condition.")

    # 4. Example Lagrangian
    print("4. Example Euclidean Lagrangian (L_E):")
    print("   For a system of non-interacting fermions, a typical Euclidean Lagrangian is:")
    print("\n   L_E = ψ̄(x,τ) * [ (∂/∂τ) + H₀ ] * ψ(x,τ)\n")
    print("   Where H₀ is the single-particle Hamiltonian (e.g., kinetic energy and chemical potential).")
    print("=" * 80)

# Execute the function to print the formula and its components
display_fermionic_partition_function()