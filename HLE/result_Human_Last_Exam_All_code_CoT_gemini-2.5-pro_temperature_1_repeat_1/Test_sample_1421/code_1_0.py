def print_fermionic_partition_function_formula():
    """
    This function prints the formula for the fermionic partition function (Z)
    in the imaginary time representation using Feynman's path integral formalism.
    """

    # Define unicode symbols for better readability
    Z = "Z"
    integral_measure = "∫ D[ψ*, ψ]"
    action = "S[ψ*, ψ]"
    beta = "β"
    tau = "τ"
    psi_star = "ψ*"
    psi = "ψ"
    H_psi = "H(ψ*, ψ)"
    d_tau = "∂_τ"
    sum_i = "Σ_i"

    # --- Print the explanation and formulas ---

    print("The formula for the fermionic partition function Z in the imaginary time path integral formalism is given by:")
    print("=" * 80)
    print(f"{Z} = {integral_measure} exp(-{action})")
    print("=" * 80)

    print("\nWhere the terms are defined as follows:")
    print(f"\n1.  {Z}: The Grand Canonical Partition Function.")
    print(f"\n2.  {integral_measure}: The path integral over all possible configurations of the fermionic")
    print(f"    (Grassmann) fields {psi} and {psi_star}.")

    print(f"\n3.  {action}: The Euclidean Action, defined as:")
    action_formula = f"    {action} = ∫[0 to {beta}] d{tau} [ {sum_i} {psi_star}_i({tau}) {d_tau} {psi}_i({tau}) + {H_psi} ]"
    print(action_formula)
    print("    - The integral is over imaginary time '{tau}' from 0 to '{beta}', where {beta} = 1/(k_B T) is the inverse temperature.".format(tau=tau, beta=beta))
    print("    - The term '{psi_star} {d_tau} {psi}' is the kinetic term.".format(psi_star=psi_star, d_tau=d_tau, psi=psi))
    print("    - '{H_psi}' is the system's Hamiltonian, where the fermion operators have been replaced by Grassmann fields.".format(H_psi=H_psi))

    print("\n4.  Boundary Condition: The path integral is subject to anti-periodic boundary conditions,")
    print("    a direct consequence of the fermionic statistics:")
    boundary_condition = f"    {psi}({beta}) = -{psi}(0)"
    print(boundary_condition)

    print("\n--- Final Combined Formula ---")
    final_formula = f"{Z} = {integral_measure} exp( -∫[0 to {beta}] d{tau} [ {sum_i} {psi_star}_i({tau}) {d_tau} {psi}_i({tau}) + {H_psi} ] )"
    print("The final equation is:")
    print(f"{Z} = {integral_measure} exp( - ∫[0 to {beta}] d{tau} [ {sum_i} {psi_star}_i({tau}) {d_tau} {psi}_i({tau}) + {H_psi} ] )")


# Execute the function to print the formula
print_fermionic_partition_function_formula()

<<<Z = ∫ D[ψ*, ψ] exp( -∫[0 to β] dτ [ Σ_i ψ*_i(τ) ∂_τ ψ_i(τ) + H(ψ*, ψ) ] )>>>