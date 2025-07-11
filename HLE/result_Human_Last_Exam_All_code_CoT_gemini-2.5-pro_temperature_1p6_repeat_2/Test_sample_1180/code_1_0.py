import sympy

def solve():
    """
    Calculates the thickness of the double point of the stable reduction.
    """
    # Define the variable t
    t = sympy.Symbol('t')

    # The stable reduction analysis leads to partitioning the branch points.
    # The curve is z^2 = 2*x^5 + 2*x^3 + 1.
    # At infinity (x=1/t), the model is w^2 = t*(t^5 + 2*t^2 + 2).
    # The branch points are the roots of t*(t^5 + 2*t^2 + 2) = 0.
    # The valuations of the roots partition the branch points into two sets.
    # The polynomials for these partitions are:
    P1 = sympy.Poly(t, t)
    P2 = sympy.Poly(t**5 + 2*t**2 + 2, t)

    # The thickness of the double point is the valuation of the resultant of these two polynomials.
    # The valuation v is the 2-adic valuation, normalized so that v(2)=1.

    # Calculate the resultant of P1 and P2
    res = sympy.resultant(P1, P2)

    # The resultant is an integer. We need to find its 2-adic valuation.
    # v(n) is the exponent of 2 in the prime factorization of n.
    resultant_val = int(res)
    
    if resultant_val == 0:
      valuation = "infinity" # Resultant is non-zero in this case
    else:
      temp_val = abs(resultant_val)
      valuation = 0
      if temp_val > 0:
        while temp_val % 2 == 0:
          valuation += 1
          temp_val //= 2
    
    print("The stable reduction is analyzed via its branch points in the t=1/x coordinate.")
    print("The branch points partition into two sets, corresponding to the polynomials:")
    print(f"P1(t) = {P1.as_expr()}")
    print(f"P2(t) = {P2.as_expr()}")
    print(f"The resultant of P1(t) and P2(t) is: {resultant_val}")
    print("The thickness of the double point is the 2-adic valuation of this resultant.")
    print("The final equation for the thickness is:")
    print(f"thickness = v({resultant_val}) = {valuation}")
    print("\nThe numbers in the final equation are:")
    print(f"Resultant: {resultant_val}")
    print(f"Thickness (valuation): {valuation}")

solve()