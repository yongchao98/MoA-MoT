def check_answer():
    """
    Checks the correctness of the answer for the given QFT problem.

    The problem asks for the mass dimension of kappa (κ) and the renormalizability
    of the theory described by the Lagrangian:
    L_int = κ * ψ_bar * σ_μν * ψ * F^μν

    The provided answer is D: [κ]_M = -1, and the theory is not renormalizable.
    """

    # Step 1: Define the known mass dimensions in 4D spacetime (natural units).
    # Dimension of the Lagrangian density L must be 4 for the action to be dimensionless.
    dim_L = 4
    # Dimension of the fermion field ψ is 3/2.
    dim_psi = 1.5
    # Dimension of the field strength tensor F^μν is 2.
    dim_F = 2
    # Dimension of σ_μν is 0, as it's composed of dimensionless gamma matrices.
    dim_sigma = 0

    # Step 2: Calculate the mass dimension of the coupling constant κ.
    # The total dimension of the interaction term must be 4.
    # The dimensions of the components in the product add up.
    # dim_L = [κ] + [ψ_bar] + [σ_μν] + [ψ] + [F^μν]
    # Note: [ψ_bar] = [ψ]
    dim_kappa = dim_L - (dim_psi + dim_sigma + dim_psi + dim_F)

    # Step 3: Determine the renormalizability based on the dimension of κ.
    # A negative mass dimension implies the theory is non-renormalizable.
    if dim_kappa < 0:
        renormalizability = "not renormalizable"
    elif dim_kappa == 0:
        renormalizability = "renormalizable"
    else: # dim_kappa > 0
        # A positive mass dimension implies super-renormalizability.
        # For the given options, this would mean the theory is renormalizable.
        renormalizability = "renormalizable"

    # Step 4: Check the calculated results against the provided answer (Option D).
    # Option D states: [κ]_M = -1, Theory is not renormalizable.
    expected_dim_kappa = -1
    expected_renormalizability = "not renormalizable"

    # Verify the mass dimension of kappa.
    if dim_kappa != expected_dim_kappa:
        return (f"Incorrect. The calculated mass dimension of kappa is {dim_kappa}, "
                f"but the answer states it is {expected_dim_kappa}.")

    # Verify the renormalizability.
    if renormalizability != expected_renormalizability:
        return (f"Incorrect. Based on the calculated dimension [κ]={dim_kappa}, "
                f"the theory is {renormalizability}, but the answer states it is "
                f"{expected_renormalizability}.")

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check and print the result.
result = check_answer()
print(result)