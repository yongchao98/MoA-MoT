import sympy

def solve_critical_exponent_order():
    """
    Symbolically determines the order in the coupling constant 'u' where the
    critical exponent nu (ν) gets its first non-trivial contribution.
    """
    # Define the symbolic variables
    # u is the coupling constant
    u = sympy.Symbol('u')
    # B is a positive constant that arises from the one-loop diagram calculation
    # for the anomalous dimension of the phi^2 operator.
    B = sympy.Symbol('B', positive=True)

    print("Step 1: Define the relationship between the critical exponent ν and the coupling u.")
    print("Within the renormalization group framework, the critical exponent ν is given by:")
    print("ν = 1 / (2 + η₂(u))")
    print("where η₂(u) is the anomalous dimension of the φ² operator.")
    print("\nPerturbative calculations (one-loop) show that the first term in η₂(u) is linear in u:")
    
    # To first order (one-loop), the anomalous dimension eta_2(u) is proportional to u.
    eta_2_function = -B * u
    print(f"η₂(u) = {eta_2_function} + O(u²)")

    # The full expression for nu as a function of u is therefore:
    nu_of_u = 1 / (2 + eta_2_function)
    print(f"\nSubstituting this into the formula for ν gives:")
    print(f"ν(u) = 1 / (2 - B*u)")

    print("\nStep 2: Expand ν(u) in a power series around u=0 to find the contributions at each order.")
    # Perform a Taylor series expansion of nu(u) around u=0.
    # We expand to the second term to see the first-order contribution.
    nu_series = nu_of_u.series(u, 0, 2)

    # Extract coefficients for clear presentation
    c0 = nu_series.coeff(u, 0)
    c1 = nu_series.coeff(u, 1)

    print(f"The Taylor series expansion of ν(u) is: ν(u) = {nu_series}")

    print("\nStep 3: Analyze the series to identify the order of the first correction.")
    print("The zeroth-order term corresponds to the mean-field value of ν, calculated at u=0.")
    print(f"ν(0) = {c0}")

    print("\nThe first correction to this value is the term proportional to u.")
    print(f"The first-order term in the expansion is: ({c1}) * u")
    print("This term is of order 1 with respect to the coupling constant u.")
    
    print("\n---")
    print("Conclusion: The initial non-vanishing contribution to ν beyond its mean-field value")
    print("is at the first order in the coupling constant u.")
    
    # As requested, printing the final equation showing each number.
    print("\nThe final equation for the expansion of ν in u is:")
    # We use pretty print for a clearer mathematical representation.
    sympy.pprint(sympy.Eq(sympy.Symbol('ν'), nu_series), use_unicode=True)

solve_critical_exponent_order()