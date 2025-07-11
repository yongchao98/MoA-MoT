import math

def P(x):
    """
    Computes the value of the polynomial P(X) at a given integer x.
    P(X) = (X^5-1)*(X^5-X)*(X^5-X^2)*(X^5-X^3)*(X^5-X^4)
    """
    if not isinstance(x, int):
        # To handle potential large numbers, ensure we're using Python's arbitrary-precision integers
        x = int(x)
    
    term1 = x**5 - 1
    term2 = x**5 - x
    term3 = x**5 - x**2
    term4 = x**5 - x**3
    term5 = x**5 - x**4
    
    # Handle the case where a term is zero to avoid issues with valuation logic
    if 0 in [term1, term2, term3, term4, term5]:
        return 0
        
    return term1 * term2 * term3 * term4 * term5

def get_valuation(n, q):
    """
    Computes the q-adic valuation of n, i.e., the exponent of the highest power of q that divides n.
    Returns a large number if n is 0, representing infinity.
    """
    if n == 0:
        return float('inf')
    if q <= 1:
        raise ValueError("q must be a prime > 1")
    
    count = 0
    while n % q == 0 and n != 0:
        count += 1
        n //= q
    return count

def find_min_valuation(q, test_range):
    """
    Finds the minimum q-adic valuation of P(x) for x in a given range,
    excluding values of x divisible by q.
    """
    min_val = float('inf')
    for x in test_range:
        if x % q == 0:
            continue
        val = get_valuation(P(x), q)
        if val < min_val:
            min_val = val
    return min_val

def main():
    """
    Main function to compute the limit L.
    """
    # Determine the exponent for q=2
    # We test odd numbers for x. A small range is sufficient to find the minimum.
    k2 = find_min_valuation(2, range(3, 20, 2))
    
    # Determine the exponent for q=3
    # We test numbers not divisible by 3.
    k3 = find_min_valuation(3, [x for x in range(2, 20) if x % 3 != 0])
    
    # Determine the exponent for q=5
    # We test numbers not divisible by 5.
    k5 = find_min_valuation(5, [x for x in range(2, 20) if x % 5 != 0])
    
    # The limit L is 2^k2 * 3^k3 * 5^k5
    val2 = 2**k2
    val3 = 3**k3
    val5 = 5**k5
    
    limit = val2 * val3 * val5
    
    print("The limit of the sequence (g_n) is L = 2^k2 * 3^k3 * 5^k5")
    print(f"The exponent for q=2 is k2 = {k2}")
    print(f"The exponent for q=3 is k3 = {k3}")
    print(f"The exponent for q=5 is k5 = {k5}")
    print(f"So the limit is {val2} * {val3} * {val5} = {limit}")

if __name__ == "__main__":
    main()