def solve():
    """
    This function calculates the required values of F(n) for the given primes.
    """

    def modInverse(n, modulus):
        """
        Calculates the modular multiplicative inverse of n modulo modulus.
        It uses the property that for a prime modulus, n^(modulus-2) is the inverse.
        """
        return pow(n, modulus - 2, modulus)

    # Case 1: p = 80039
    p1 = 80039
    # The final expression simplifies to S(3) = 7/8
    num1 = 7
    den1 = 8
    
    # Calculate 7/8 mod 80039
    den1_inv = modInverse(den1, p1)
    ans1 = (num1 * den1_inv) % p1

    # Case 2: p = 80077
    p2 = 80077
    # The final expression simplifies to 27/4
    num2 = 27
    den2 = 4

    # Calculate 27/4 mod 80077
    den2_inv = modInverse(den2, p2)
    ans2 = (num2 * den2_inv) % p2

    print(f"For p = {p1}, the final equation is ({num1} * ({den1}^-1)) mod {p1}")
    print(f"The result is ({num1} * {den1_inv}) mod {p1} = {ans1}")
    print("-" * 20)
    print(f"For p = {p2}, the final equation is ({num2} * ({den2}^-1)) mod {p2}")
    print(f"The result is ({num2} * {den2_inv}) mod {p2} = {ans2}")
    print("-" * 20)
    print(f"The two answers separated by a comma are: {ans1},{ans2}")

solve()
<<<70035,20020>>>