import sympy

def solve_critical_exponent():
    """
    Calculates the critical exponent p for the L^p integrability of the integral I
    using a scaling argument.
    """

    # Define the variables of integration
    x, y = sympy.symbols('x y')

    # List of monomials in the phase of the integral
    monomials = [
        x, y,
        x**2, x*y, y**2,
        x**3, x**2*y, x*y**2, y**3
    ]

    # Dimension of the integration domain [0,1]^2
    n = 2
    print(f"The dimension of the integration domain (dxdy) is n = {n}.")
    print("-" * 30)

    # Calculate the degree of each monomial
    degrees = [sympy.total_degree(m, x, y) for m in monomials]
    
    print("The monomials in the exponent correspond to the parameters a_1, ..., a_9.")
    print("Their respective degrees are:")
    print(degrees)
    print("-" * 30)

    # Sum of the degrees
    D = sum(degrees)
    print(f"The sum of the degrees is D = {D}.")
    print("-" * 30)
    
    print("The scaling argument leads to the equation relating p, n, and D:")
    print(f"-n * p + D = 0")
    print("Substituting the values of n and D, we get:")
    print(f"-{n} * p + {D} = 0")
    print("-" * 30)

    # Solve for the critical exponent p
    # The equation is -n*p + D = 0
    p = D / n

    print("Solving for p:")
    print(f"{n} * p = {D}")
    print(f"p = {D} / {n}")
    print(f"p = {p}")
    print("-" * 30)
    
    print(f"The critical exponent is {p}.")
    print(f"The function I is in L^q for q > {p} and not in L^q for q <= {p}.")
    print(f"Thus, the largest p for which I is not in L^p(R^9) is {p}.")

solve_critical_exponent()