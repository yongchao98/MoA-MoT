import math

def solve_integral():
    """
    Calculates the value of the given mathematical expression based on analytical derivation.
    """
    pi = math.pi

    # The integral was analytically solved and simplified to the expression:
    # (8/15)*pi^8 + (1/3)*pi^2 - (1/2)*pi + 1
    # Here, we calculate each term of this final expression.

    # Term 1: from the p^7 component of the integral
    term1 = (8 / 15) * (pi**8)

    # Term 2: from the p and pe^-p components
    term2 = (1 / 3) * (pi**2)

    # Term 3: from the sinh component
    term3 = (-1 / 2) * pi
    
    # Term 4: the constant component
    term4 = 1

    print("The final expression is a sum of four terms:")
    print(f"Equation: (8/15)*pi^8 + (1/3)*pi^2 - (1/2)*pi + 1")
    print("-" * 30)
    print(f"Term 1 value: {term1}")
    print(f"Term 2 value: {term2}")
    print(f"Term 3 value: {term3}")
    print(f"Term 4 value: {term4}")
    print("-" * 30)
    
    # Calculate the final numerical value
    total_value = term1 + term2 + term3 + term4
    
    print(f"Final calculated value: {total_value}")

solve_integral()