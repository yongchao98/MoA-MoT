def display_fermionic_partition_function_formula():
    """
    Prints the formula for the fermionic partition function Z
    in the imaginary time path integral formalism.
    """

    print("The formula for the fermionic partition function (Z) in the imaginary time representation using Feynman’s path integral formalism is as follows:")
    print("-" * 120)

    # Main formula for Z
    z_formula = "Z = ∫ Dψ Dψ̄  exp(-S[ψ, ψ̄])"
    print("\n[1] The Partition Function Formula:")
    print(f"    {z_formula}")

    # Definition of the Euclidean Action S
    action_formula = "S[ψ, ψ̄] = ∫_0^β dτ [ ψ̄(τ) (∂_τ) ψ(τ) + H(ψ̄(τ), ψ(τ)) ]"
    print("\n[2] The Euclidean Action (S):")
    print(f"    {action_formula}")

    # Boundary Conditions
    boundary_condition_intro = "The path integral is performed over all fields satisfying anti-periodic boundary conditions:"
    psi_bc = "ψ(β) = -ψ(0)"
    print("\n[3] The Boundary Conditions:")
    print(f"    {boundary_condition_intro}")
    print(f"    {psi_bc}")

    # Highlighting numbers in the equations as requested
    print("\n" + "-" * 120)
    print("Analysis of the components and numbers in the formula:")
    print("\nIn the Action formula [2]:")
    print("  - The integral ∫_0^β has a lower limit, which is the number 0.")
    print("  - The upper limit is β, which represents the inverse temperature (1/k_B T).")

    print("\nIn the Boundary Condition formula [3]:")
    print("  - The equation ψ(β) = -ψ(0) explicitly defines the state at time τ = 0 and τ = β.")
    print("  - The negative sign '-' corresponds to the number -1, which is characteristic of fermions.")
    print("    It can be written as: ψ(β) = (-1) * ψ(0).")

if __name__ == '__main__':
    display_fermionic_partition_function_formula()