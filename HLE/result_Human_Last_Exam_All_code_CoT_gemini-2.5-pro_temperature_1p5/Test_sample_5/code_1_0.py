def solve_gamma_product():
    """
    Calculates the proportionality factor C in the equation:
    gamma_{mu nu} gamma_{mu_1 ... mu_k} gamma^{mu nu} = C * gamma_{mu_1 ... mu_k}
    """
    
    # Let d be the number of spacetime dimensions and k be the rank of the gamma matrix tensor.
    # These are symbolic variables.
    d_var = "d"
    k_var = "k"

    print("The problem is to find the proportionality factor C in the equation:")
    print(f"γ_{{μν}} γ_{{μ₁…μₖ}} γ^{{μν}} = C γ_{{μ₁…μₖ}}\n")

    print("We use the general identity for products of gamma matrices:")
    print(f"γ^{{ρ₁…ρₚ}} γ_{{μ₁…μₖ}} γ_{{ρ₁…ρₚ}} = (-1)^{{pk + p(p-1)/2}} p! (d-k)!/(p!(d-k-p)!) γ_{{μ₁…μₖ}}\n")
    
    print("In our case, the outer gamma matrices form a rank-2 tensor, so p = 2.")
    
    # Step 1: Calculate the sign factor
    sign_factor_expr = "(-1)^{{2k + 2(2-1)/2}}"
    sign_factor_simpl = "(-1)^{{2k + 1}}"
    sign_factor_final = -1
    print(f"1. The sign factor is {sign_factor_expr} = {sign_factor_simpl} = {sign_factor_final}.")
    
    # Step 2: Calculate the factorial term
    p_factorial_val = 2
    print(f"2. The factorial term p! for p=2 is 2! = {p_factorial_val}.")

    # Step 3: The binomial coefficient
    binom_coeff_expr = f"C({d_var}-{k_var}, 2)"
    binom_coeff_expanded = f"({d_var}-{k_var})({d_var}-{k_var}-1) / 2"
    print(f"3. The binomial coefficient is C(d-k, p) = {binom_coeff_expr} = {binom_coeff_expanded}.")

    # Step 4: Combine the terms to find the proportionality factor C
    print("\nMultiplying these components together gives the proportionality factor C:")
    final_C_expr = f"C = ({sign_factor_final}) * ({p_factorial_val}) * ({binom_coeff_expanded})"
    final_C_simpl = f"-({d_var}-{k_var})({d_var}-{k_var}-1)"
    
    print(final_C_expr)
    print(f"C = {final_C_simpl}\n")

    print("Therefore, the final relation is:")
    
    term1 = f"({d_var}-{k_var})"
    term2 = f"({d_var}-{k_var}-1)"
    
    print(f"γ_{{μν}} γ_{{μ₁…μₖ}} γ^{{μν}} = -{term1}{term2} γ_{{μ₁…μₖ}}")


solve_gamma_product()