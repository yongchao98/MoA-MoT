import sympy

def explain_critical_exponent_order():
    """
    Explains and calculates the order of the first correction to the critical exponent nu.
    """

    # --- Step 1: Theoretical Background ---
    print("In perturbative renormalization group analysis of phi^4 theory, critical exponents are computed as an expansion in the coupling constant 'u'.")
    print("We want to find the order of the first non-vanishing correction to the mean-field value of the critical exponent nu (ν).\n")

    # --- Step 2: The Expansion of nu(u) ---
    print("The exponent ν can be expressed as a power series in the coupling u around the Gaussian fixed point (u=0):")
    print("ν(u) = ν(0) + c_1*u^1 + c_2*u^2 + ...")
    print("The mean-field value is ν(0) = 1/2.")
    print("The first non-vanishing correction is the c_1*u^1 term, as the one-loop calculation yields a non-zero coefficient c_1.\n")
    
    # --- Step 3: Derivation in the epsilon-expansion ---
    print("This is explicitly shown in the ε-expansion (ε = 4-d).")
    print("The non-trivial fixed point coupling u* is found to be of order ε, i.e., u* ∝ ε.")
    print("The formula for ν for the O(n)-symmetric model, calculated to the first order in ε, is:")
    print("ν = 1 / (2 - (n+2)/(n+8) * ε)\n")

    # --- Step 4: Show the specific calculation for n=1 (Ising Model) ---
    n = 1
    nu_mf = sympy.Rational(1, 2)
    coeff_num = n + 2
    coeff_den = n + 8
    
    print("Let's look at the specific case for the Ising-like model where n = 1.")
    print(f"ν = 1 / (2 - ({n}+{2})/({n}+{8}) * ε)")
    print(f"ν = 1 / (2 - {coeff_num}/{coeff_den} * ε)")
    
    print("\nExpanding this for small ε gives:")
    print(f"ν ≈ (1/2) * (1 + (1/2) * ({coeff_num}/{coeff_den}) * ε)")
    
    final_coeff_num = coeff_num
    final_coeff_den = 2 * coeff_den
    
    # Use sympy to simplify the fraction
    final_coeff_simplified = sympy.Rational(final_coeff_num, final_coeff_den)

    print(f"ν ≈ {nu_mf} + ({final_coeff_num} / {final_coeff_den}) * ε")
    print(f"ν ≈ {nu_mf} + {final_coeff_simplified} * ε\n")
    
    # --- Step 5: Conclusion ---
    order = 1
    print(f"Since the first correction term is proportional to ε, and u* is proportional to ε, the correction is linear in the coupling u*.")
    print(f"Therefore, the initial non-vanishing contribution to ν appears at order {order} in the coupling constant u.\n")

# Run the explanation
explain_critical_exponent_order()
