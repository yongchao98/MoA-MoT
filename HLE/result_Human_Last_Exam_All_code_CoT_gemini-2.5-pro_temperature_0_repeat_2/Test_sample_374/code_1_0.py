def solve():
    """
    Calculates the highest possible order for the inertial quotient E.
    This is equivalent to calculating the order of the group GL(4, 2).
    """
    n = 4
    q = 2

    # The order of GL(n, q) is the product of (q^n - q^i) for i from 0 to n-1.
    
    # Calculate the terms in the product
    term1 = q**n - q**0
    term2 = q**n - q**1
    term3 = q**n - q**2
    term4 = q**n - q**3
    
    # Calculate the final order
    order = term1 * term2 * term3 * term4

    print("The highest possible order for the inertial quotient E is the order of GL(4, 2).")
    print("The calculation is based on the formula for |GL(n, q)| with n=4 and q=2.")
    print("\nThe terms of the product are:")
    print(f"(2^4 - 2^0) = (16 - 1) = {term1}")
    print(f"(2^4 - 2^1) = (16 - 2) = {term2}")
    print(f"(2^4 - 2^2) = (16 - 4) = {term3}")
    print(f"(2^4 - 2^3) = (16 - 8) = {term4}")
    
    print("\nThe final equation is:")
    print(f"{term1} * {term2} * {term3} * {term4} = {order}")
    
    print(f"\nThus, the highest order that E can have is {order}.")

solve()