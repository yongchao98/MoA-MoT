import math

def solve():
    """
    Calculates the highest possible order for the inertial quotient E.

    This corresponds to the order of the outer automorphism group of the
    defect group D, which is an elementary abelian group of order 16.
    """

    # The defect group D is elementary abelian of order 16.
    # This means D is a 4-dimensional vector space over the field F_2.
    # The automorphism group Aut(D) is GL(4, F_2).
    # Since D is abelian, Out(D) = Aut(D).
    # The highest order for E is |GL(4, F_2)|.

    n = 4  # Dimension of the vector space
    q = 2  # Size of the finite field

    # The order of GL(n, q) is product_{i=0 to n-1} (q^n - q^i)
    
    term1 = q**n - q**0
    term2 = q**n - q**1
    term3 = q**n - q**2
    term4 = q**n - q**3

    # Calculate the product of the terms
    order = term1 * term2 * term3 * term4
    
    print("The highest possible order for E is the order of GL(4, F_2).")
    print("The calculation for the order is:")
    print(f"({q**n} - {q**0}) * ({q**n} - {q**1}) * ({q**n} - {q**2}) * ({q**n} - {q**3})")
    print(f"= {term1} * {term2} * {term3} * {term4}")
    print(f"= {order}")

solve()
