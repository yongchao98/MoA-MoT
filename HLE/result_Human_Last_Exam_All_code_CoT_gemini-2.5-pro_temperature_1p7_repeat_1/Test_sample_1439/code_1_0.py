import sympy

def solve_critical_exponent_order():
    """
    Determines the order in the coupling constant u at which the critical
    exponent ν first receives a non-vanishing contribution in φ⁴ theory.
    """

    # Define the symbolic variables
    # u: the coupling constant
    # c1, c2: coefficients from one-loop and two-loop calculations
    u = sympy.Symbol('u')
    c1 = sympy.Symbol('c1')
    c2 = sympy.Symbol('c2')

    print("Step 1: The relationship between ν and γ₂(u)")
    print("-------------------------------------------------")
    print("In the renormalization group analysis of φ⁴ theory, the critical exponent ν")
    print("is given by the formula:")
    print("  ν(u) = 1 / (2 - γ₂(u))")
    print("where u is the coupling constant and γ₂(u) is the anomalous dimension of")
    print("the φ² operator.\n")

    print("Step 2: The perturbative expansion of γ₂(u)")
    print("-------------------------------------------------")
    print("The anomalous dimension γ₂(u) is calculated as a power series in u. The")
    print("leading contribution comes from a one-loop diagram, which is proportional to u.")
    print("So, the series expansion starts as:")
    gamma_2_expansion = c1*u + c2*u**2 + sympy.O(u**3)
    print(f"  γ₂(u) = {sympy.printing.pretty(gamma_2_expansion, use_unicode=False)}\n")
    print("(Note: c1 is non-zero.)\n")


    print("Step 3: Taylor expansion of ν(u)")
    print("-------------------------------------------------")
    print("To find the first correction to the mean-field value (ν at u=0), we")
    print("expand ν(u) in a power series around u = 0.")
    
    # The expression for nu(u)
    nu_expression = 1 / (2 - gamma_2_expansion)

    # Perform the series expansion
    nu_series = nu_expression.series(u, 0, 3)

    print("The series expansion for ν(u) is:")
    # We will print the equation clearly showing the numbers
    # The sympy output for nu_series is: 1/2 + c1*u/4 + (c1**2/16 + c2/8)*u**2 + O(u**3)
    print("  ν(u) = 1/2 + (c1/4)*u + ...\n")
    

    print("Step 4: Conclusion")
    print("-------------------------------------------------")
    print("The expansion ν(u) = 1/2 + (c1/4)*u + O(u²) shows that:")
    print(" - The term of order u⁰ is 1/2, which is the classical mean-field value.")
    print(" - The first non-vanishing contribution (correction) to this value is the term")
    print("   proportional to u¹ (the first power of u).")
    print("\nTherefore, the critical exponent ν acquires its initial non-vanishing contribution")
    print("at the first order in the coupling constant u.")

solve_critical_exponent_order()
<<<1>>>