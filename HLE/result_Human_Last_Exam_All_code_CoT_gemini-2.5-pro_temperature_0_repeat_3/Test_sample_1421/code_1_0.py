def display_fermionic_partition_function():
    """
    Prints the formula for the fermionic partition function Z
    in the imaginary time path integral formalism.
    """

    # The main formula for the partition function Z
    partition_function_formula = (
        "Z = ∫ Dψ̄Dψ exp(-S_E[ψ̄, ψ])"
    )

    # The formula for the Euclidean action S_E
    # Here, H is the Hamiltonian operator for the system.
    # The integral is over imaginary time τ from 0 to β, and over d spatial dimensions x.
    euclidean_action_formula = (
        "S_E[ψ̄, ψ] = ∫[from τ=0 to β] dτ ∫ d^d x  ψ̄(x, τ) * (∂/∂τ + H) * ψ(x, τ)"
    )

    # The anti-periodic boundary condition for fermionic fields
    # This is a consequence of the trace over the Fock space of fermions.
    boundary_condition = (
        "Boundary Condition: ψ(x, τ=0) = -ψ(x, τ=β)"
    )

    print("The formula for the fermionic partition function (Z) in the imaginary time path integral formalism is given by:")
    print("-" * 80)
    print("1. Partition Function (Z):")
    print(f"   {partition_function_formula}\n")

    print("2. Where S_E is the Euclidean Action:")
    print(f"   {euclidean_action_formula}\n")

    print("3. And the Grassmann fields (ψ, ψ̄) obey the anti-periodic boundary condition:")
    print(f"   {boundary_condition}")
    print("-" * 80)
    print("\nHere,")
    print("  - Z: The partition function.")
    print("  - ∫ Dψ̄Dψ: The functional (or path) integral over all configurations of the anti-commuting Grassmann fields ψ and ψ̄.")
    print("  - S_E: The Euclidean action of the system.")
    print("  - β: Inverse temperature (1 / (k_B * T)).")
    print("  - H: The Hamiltonian of the system.")
    print("  - The numbers in the formula are the integral limits for imaginary time τ, which are 0 and β.")


if __name__ == "__main__":
    display_fermionic_partition_function()