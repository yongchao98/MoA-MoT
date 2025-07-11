import math

def solve():
    """
    Calculates the limit of the sequence g_n based on the polynomial P(X).

    The limit g is the greatest common divisor of P(p) for all sufficiently large primes p.
    Our analysis showed that the prime factors of g can only be 2, 3, and 5.
    We then determined the exponents of these prime factors in g by finding the minimum
    p-adic valuation of P(p) over all large primes p.

    The exponents found are:
    - For prime 2: The minimum valuation v_2(P(p)) is 10.
    - For prime 3: The minimum valuation v_3(P(p)) is 2.
    - For prime 5: The minimum valuation v_5(P(p)) is 1.

    Thus, the limit g = 2^10 * 3^2 * 5^1.
    """
    
    # The exponents determined from the mathematical analysis.
    exp_2 = 10
    exp_3 = 2
    exp_5 = 1
    
    # The prime bases
    base_2 = 2
    base_3 = 3
    base_5 = 5
    
    # Calculate the components
    term1 = base_2 ** exp_2
    term2 = base_3 ** exp_3
    term3 = base_5 ** exp_5
    
    # Calculate the final limit g
    g = term1 * term2 * term3
    
    # Print the explanation and the result as an equation
    print("Based on the analysis, the limit g is calculated as the product of prime powers:")
    print(f"g = {base_2}^{exp_2} * {base_3}^{exp_3} * {base_5}^{exp_5}")
    print(f"g = {term1} * {term2} * {term3}")
    print(f"The limit of the sequence is {g}")

solve()