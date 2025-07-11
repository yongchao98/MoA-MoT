def solve():
    """
    Solves the problem for the two given primes.
    """

    # Helper function for modular inverse
    # In Python 3.8+, pow(a, -1, m) can be used.
    def mod_inverse(a, m):
        return pow(a, m - 2, m)

    # Case 1: p = 80039
    p1 = 80039
    print(f"--- Calculating for p = {p1} ---")
    
    # As derived in the explanation, the problem reduces to calculating S(3) mod p.
    # S(3) = 7/8
    num1 = 7
    den1 = 8
    print(f"The problem reduces to calculating the value of S(n) for n=3, which is {num1}/{den1}.")
    
    # Calculate 7 * (8^-1) mod 80039
    inv_den1 = mod_inverse(den1, p1)
    result1 = (num1 * inv_den1) % p1
    print(f"The value is ({num1} * {den1}^-1) mod {p1} = {result1}")
    print("")

    # Case 2: p = 80077
    p2 = 80077
    print(f"--- Calculating for p = {p2} ---")
    
    # As derived, the problem reduces to calculating S(5-2p) mod p.
    # S(5-2p) evaluates to 27/4.
    num2 = 27
    den2 = 4
    print(f"The problem reduces to calculating a value that simplifies to {num2}/{den2}.")

    # Calculate 27 * (4^-1) mod 80077
    inv_den2 = mod_inverse(den2, p2)
    result2 = (num2 * inv_den2) % p2
    print(f"The value is ({num2} * {den2}^-1) mod {p2} = {result2}")
    print("")

    # Final Answer
    print(f"The final answers are separated by a comma:")
    print(f"{result1},{result2}")

solve()
<<<70035,20026>>>