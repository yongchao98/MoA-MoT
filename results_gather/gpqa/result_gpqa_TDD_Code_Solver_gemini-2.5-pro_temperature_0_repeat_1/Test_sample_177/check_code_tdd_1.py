def check_theory_properties(option: str):
    """
    Checks the mass dimension of kappa and the renormalizability of the theory
    for the given Lagrangian L_int = kappa * bar(psi) * sigma_mu_nu * psi * F^mu_nu.

    Args:
        option (str): The selected option ('A', 'B', 'C', or 'D').

    Returns:
        str: "Correct" if the option is correct, otherwise a reason for the error.
    """
    # --- Step 1: Define known mass dimensions in natural units (hbar=c=1) ---
    # The action S = integral(d^4x * L) must be dimensionless ([S]=0).
    # The dimension of the spacetime volume element d^4x is [d^4x] = -4.
    # Therefore, the Lagrangian density L must have a mass dimension of [L] = 4.
    dim_L = 4

    # The kinetic term for a fermion is L_fermion = bar(psi) * (i * gamma^mu * partial_mu - m) * psi.
    # The dimension of the derivative term must match [L]=4.
    # [bar(psi)] + [partial_mu] + [psi] = 4
    # Since [bar(psi)] = [psi] and [partial_mu] = 1, we have:
    # 2 * [psi] + 1 = 4  =>  [psi] = 3/2
    dim_psi = 1.5

    # The kinetic term for the electromagnetic field is L_EM = -1/4 * F_munu * F^munu.
    # The dimension must be [L]=4.
    # [F_munu] + [F^munu] = 4  =>  2 * [F_munu] = 4  =>  [F_munu] = 2
    dim_F_munu = 2

    # The term sigma_munu = (i/2) * [gamma_mu, gamma_nu] is composed of dimensionless
    # gamma matrices, so it is also dimensionless.
    dim_sigma_munu = 0

    # --- Step 2: Calculate the mass dimension of the coupling constant kappa ---
    # The interaction term L_int = kappa * bar(psi) * sigma_munu * psi * F^munu must also have [L_int] = 4.
    # [L_int] = [kappa] + [bar(psi)] + [sigma_munu] + [psi] + [F^munu]
    # 4 = [kappa] + [psi] + [sigma_munu] + [psi] + [F_munu]
    # 4 = [kappa] + 1.5 + 0 + 1.5 + 2
    # 4 = [kappa] + 5
    calculated_dim_kappa = dim_L - (dim_psi + dim_sigma_munu + dim_psi + dim_F_munu)
    # calculated_dim_kappa = 4 - (1.5 + 0 + 1.5 + 2) = 4 - 5 = -1

    # --- Step 3: Determine renormalizability ---
    # A theory is non-renormalizable by power-counting if its coupling constant
    # has a negative mass dimension.
    is_renormalizable_calculated = (calculated_dim_kappa >= 0)
    renormalizability_calculated_str = "renormalizable" if is_renormalizable_calculated else "not renormalizable"

    # --- Step 4: Define the claims of each option ---
    options = {
        'A': {'dim': -1, 'renorm': True, 'renorm_str': 'renormalizable'},
        'B': {'dim': 1, 'renorm': True, 'renorm_str': 'renormalizable'},
        'C': {'dim': 1, 'renorm': False, 'renorm_str': 'not renormalizable'},
        'D': {'dim': -1, 'renorm': False, 'renorm_str': 'not renormalizable'}
    }

    if option.upper() not in options:
        return f"Invalid option '{option}'. Please choose from A, B, C, or D."

    selected_option_claims = options[option.upper()]
    claimed_dim = selected_option_claims['dim']
    claimed_renorm = selected_option_claims['renorm']
    claimed_renorm_str = selected_option_claims['renorm_str']

    # --- Step 5: Compare calculated values with the claims of the selected option ---
    # Check mass dimension
    if claimed_dim != calculated_dim_kappa:
        return (f"Incorrect: The mass dimension of kappa is incorrect. "
                f"Calculation: [L_int] = [kappa] + 2*[psi] + [F_munu] => 4 = [kappa] + 2*(3/2) + 2 = [kappa] + 5. "
                f"Therefore, the calculated [kappa] = -1. The answer claims [kappa] = {claimed_dim}.")

    # Check renormalizability
    if claimed_renorm != is_renormalizable_calculated:
        return (f"Incorrect: The conclusion on renormalizability is wrong. "
                f"Since the mass dimension of the coupling constant [kappa] = {calculated_dim_kappa} is negative, "
                f"the theory is non-renormalizable by power counting. The answer claims the theory is {claimed_renorm_str}.")

    # If both checks pass, the option is correct.
    return "Correct"

# Example of how to use the checker for each option:
# The provided answer is a placeholder, so we can't check it directly.
# Instead, we can determine the correct option by testing all possibilities.
print(f"Checking option A: {check_theory_properties('A')}")
print(f"Checking option B: {check_theory_properties('B')}")
print(f"Checking option C: {check_theory_properties('C')}")
print(f"Checking option D: {check_theory_properties('D')}")