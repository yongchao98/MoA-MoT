def display_fermionic_partition_function_formula():
    """
    This function prints the formula for the fermionic partition function
    in the imaginary time path integral formalism.
    """

    # --- Define symbols for the formula ---
    partition_function = "Z"
    integral_measure = "∫ Dψ̄ Dψ"
    exponential = "exp"
    euclidean_action = "S_E[ψ̄, ψ]"

    # --- Print the main formula for the partition function Z ---
    print("The fermionic partition function Z is given by:")
    formula_z = f"{partition_function} = {integral_measure} {exponential}(-{euclidean_action})"
    print(formula_z)
    print("-" * 30)

    # --- Define and print the Euclidean action S_E ---
    print("Where the Euclidean action S_E is:")
    integral_time = "∫_0^β dτ"  # Integration from 0 to beta
    integral_space = "∫ dᵈx"
    lagrangian_euclidean = "L_E(ψ̄(x,τ), ψ(x,τ))"
    formula_s = f"{euclidean_action} = {integral_time} {integral_space} {lagrangian_euclidean}"
    print(formula_s)
    print("-" * 30)

    # --- State the crucial boundary conditions for fermions ---
    print("The path integral is taken over all Grassmann fields ψ and ψ̄")
    print("subject to anti-periodic boundary conditions in imaginary time:")
    boundary_condition = "ψ(x, 0) = -ψ(x, β)"
    print(boundary_condition)

if __name__ == '__main__':
    display_fermionic_partition_function_formula()

<<<Z = ∫ Dψ̄ Dψ exp(-S_E[ψ̄, ψ]), where S_E[ψ̄, ψ] = ∫_0^β dτ ∫ dᵈx L_E(ψ̄(x,τ), ψ(x,τ)) with anti-periodic boundary conditions ψ(x, 0) = -ψ(x, β)>>>