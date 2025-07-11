def solve_involutions():
    """
    Calculates the number of involutions in PSU(4, 997).
    """
    q = 997

    # The formula for the number of involutions in PSU(4, q) for q = 1 mod 4 is:
    # N = (q^5 * (q+1)^2 * (q^3 + 1)) / 2
    
    # Let's calculate each part of the formula.
    # Python's integers handle arbitrarily large numbers, so overflow is not an issue.
    term1 = q**5
    term2 = (q + 1)**2
    term3 = q**3 + 1
    
    # Calculate the numerator
    numerator = term1 * term2 * term3
    
    # The result is an integer, so we use integer division.
    num_involutions = numerator // 2
    
    print("This script calculates the number of involutions in PSU(4, q) where q = 997.")
    print("The formula used is N = (q^5 * (q+1)^2 * (q^3+1)) / 2.")
    print(f"For q = {q}, the components of the numerator are:")
    print(f"q^5 = {term1}")
    print(f"(q+1)^2 = {term2}")
    print(f"q^3+1 = {term3}")
    
    print("\nThe final equation with these numbers is:")
    print(f"({term1} * {term2} * {term3}) / 2 = {num_involutions}")

solve_involutions()
