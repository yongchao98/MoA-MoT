import sympy

def main():
    """
    This script explains the derivation for finding the order of the first
    correction to the critical exponent ν in φ⁴ theory.
    """
    # Define symbols for the explanation
    u = sympy.Symbol('u')
    N = sympy.Symbol('N')
    
    print("In the ε-expansion framework for φ⁴ theory, we analyze the critical exponents as a power series in the coupling constant u.")
    print("-" * 50)

    # Step 1: Relation for ν
    print("Step 1: The critical exponent ν is related to the anomalous dimension of the φ² operator, η_φ²(u).")
    print("The formula is:")
    print("1/ν = 2 - η_φ²(u)")
    print("\nIn the free theory (u=0), η_φ²(u) is zero. This gives the mean-field value:")
    print("1/ν = 2 - 0  =>  ν = 1/2")
    print("-" * 50)

    # Step 2: Expansion of the anomalous dimension
    print("Step 2: The anomalous dimension η_φ²(u) is calculated perturbatively.")
    print("The one-loop diagram gives the first term in its expansion. For an N-component theory,")
    print("in a common renormalization scheme, this is:")
    print("η_φ²(u) = (N+2)/2 * u + O(u²)")
    print("\nThe crucial fact is that this expansion starts at the first order in u (u¹).")
    print("-" * 50)
    
    # Step 3: Expanding ν in powers of u
    print("Step 3: Substitute the expansion of η_φ²(u) into the formula for ν.")
    # In a common scheme, the constant C for the first term is (N+2)/2. We'll use a generic 'C'.
    C = sympy.Symbol('C')
    eta_phi_sq = C * u + sympy.O(u**2)
    nu_inv_expr = 2 - eta_phi_sq
    
    print(f"1/ν(u) = {nu_inv_expr}")
    
    nu_expr = 1 / (2 - C*u)
    nu_series = nu_expr.series(u, 0, 3)
    
    print("\nTo find ν(u), we invert this and perform a Taylor expansion around u=0:")
    print(f"ν(u) = (2 - C*u - ... )⁻¹")
    print(f"ν(u) ≈ {nu_series}")
    print("-" * 50)

    # Step 4: Conclusion
    print("Step 4: Analyze the resulting series for ν(u).")
    print(f"The expansion ν(u) = 1/2 + (C/4)*u + ... shows:")
    print(" - A constant term (order u⁰), which is the mean-field value: 1/2")
    print(" - The first correction term, which is proportional to: u¹")
    
    final_order = 1
    mean_field_numerator = 1
    mean_field_denominator = 2

    print("\nTherefore, the initial non-vanishing contribution to ν (beyond its mean-field value) occurs at the first order in the coupling constant u.")
    print("\nFinal Equation Breakdown:")
    print(f"The equation for ν starts as: ν(u) = {mean_field_numerator}/{mean_field_denominator} + (constant) * u^{final_order} + ...")
    print(f"The numbers in this conceptual equation are:")
    print(f"Numerator of mean-field value: {mean_field_numerator}")
    print(f"Denominator of mean-field value: {mean_field_denominator}")
    print(f"Order of the first correction: {final_order}")

if __name__ == "__main__":
    main()