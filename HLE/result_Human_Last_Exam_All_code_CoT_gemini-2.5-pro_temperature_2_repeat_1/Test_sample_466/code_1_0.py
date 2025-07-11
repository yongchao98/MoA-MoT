import math

def solve():
    """
    Calculates the number of internal adjunctions in the given 2-category from F_11^3 to itself.
    This number is equivalent to the order of the general linear group GL(3, F_11).
    """
    # Parameters for the finite field and vector space dimension
    q = 11
    n = 3

    # Calculate terms for the order of GL(n, F_q)
    # The formula is |GL(n, F_q)| = (q^n - q^0) * (q^n - q^1) * ... * (q^n - q^(n-1))
    
    # Calculate q^n once to save computation
    q_to_n = q**n

    # Calculate each term in the product
    term1_val = q_to_n - (q**0)
    term2_val = q_to_n - (q**1)
    term3_val = q_to_n - (q**2)
    
    # Calculate the final result
    result = term1_val * term2_val * term3_val
    
    print(f"The number of internal adjunctions is the order of the general linear group GL({n}, F_{q}).")
    print(f"The calculation formula is: (|GL({n}, F_{q})| = ({q}^{n} - {q}^0) * ({q}^{n} - {q}^1) * ({q}^{n} - {q}^2))")
    print(f"Plugging in the values, we get: ({q_to_n} - {q**0}) * ({q_to_n} - {q**1}) * ({q_to_n} - {q**2})")
    print(f"This simplifies to the product: {term1_val} * {term2_val} * {term3_val}")
    # Using math.trunc to get a clean integer representation
    print(f"The final number of internal adjunctions is: {math.trunc(result)}")

solve()