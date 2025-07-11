import math

def solve():
    """
    Solves the problem based on the provided values of p, k, and n.
    """
    try:
        p_str = input("Enter the odd prime p (>=3): ")
        p = int(p_str)
        if p < 3 or p % 2 == 0:
            print("p must be an odd prime >= 3.")
            return
        
        # A simple primality test for demonstration.
        is_prime = True
        if p > 3:
            for i in range(2, int(math.sqrt(p)) + 1):
                if p % i == 0:
                    is_prime = False
                    break
        if not is_prime:
            print(f"{p} is not a prime number.")
            return

        k_str = input("Enter the integer k (>=1): ")
        k = int(k_str)
        if k < 1:
            print("k must be an integer >= 1.")
            return

        n_str = input("Enter the integer n (>=1): ")
        n = int(n_str)
        if n < 1:
            print("n must be an integer >= 1.")
            return
    except ValueError:
        print("Invalid input. Please enter integers.")
        return

    # The leading coefficient of the final polynomial is 1.
    # Since p >= 3 and k >= 1, p^k >= 3. So 1 is never divisible by p^k.
    # This guarantees at least one coefficient is not divisible by p^k.
    num_coeffs = 1
    
    # We analyze the coefficient of x^1, c_1. Its p-adic valuation is:
    # v_p(c_1) = p^n * ( (p^k - 1)/(p-1) - k )
    
    vp_term = (p**k - 1) // (p - 1) - k
    vp_c1 = p**n * vp_term
    
    # We check if c_1 is divisible by p^k, which is equivalent to checking if v_p(c_1) >= k.
    if vp_c1 < k:
        # This case only happens for k=1.
        num_coeffs += 1
    
    print(f"\nFor p={p}, k={k}, n={n}:")
    print(f"The number of coefficients not divisible by p^k is {num_coeffs}.")


solve()