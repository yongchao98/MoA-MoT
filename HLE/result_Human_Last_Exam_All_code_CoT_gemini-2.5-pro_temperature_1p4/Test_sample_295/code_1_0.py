def compute_reduced_euler_characteristic_mod_k(k):
    """
    Computes the reduced Euler characteristic of a specific simplicial complex,
    Delta_k, modulo a prime k >= 3.

    The problem states that k is a prime k >= 3.
    This implies k is an odd prime.

    Args:
        k: An odd prime number.
    """
    if not (isinstance(k, int) and k >= 3):
        print(f"Error: k must be an integer prime >= 3. Received k={k}.")
        return

    # The derived formula for the result is (3-k)/2 mod k.
    # Since k is an odd prime, k is odd.
    # Therefore, 3-k is an even number, and the division yields an integer.
    numerator = 3 - k
    value = numerator // 2

    # The result is the value modulo k.
    # Python's % operator for a % n returns a result with the same sign as the
    # divisor n, which works perfectly for finding the canonical representative in Z_k.
    result_mod_k = value % k

    print(f"For the prime k = {k}:")
    # Using Unicode characters for mathematical notation
    # χ̂(Δₖ) ≡ (3 - k) / 2 (mod k)
    print(f"The formula is: \u03C7\u0302(\u0394_{k}) \u2261 (3 - k) / 2 (mod {k})")
    print(f"Step 1: Calculate the numerator -> 3 - {k} = {numerator}")
    print(f"Step 2: Perform the integer division -> {numerator} // 2 = {value}")
    print(f"Step 3: Compute the modulus -> {value} % {k} = {result_mod_k}")
    print(f"The final result for k={k} is: {result_mod_k}")
    print("-" * 30)

if __name__ == "__main__":
    # A list of small prime numbers to demonstrate the calculation
    primes = [3, 5, 7, 11, 13, 17, 19]
    for p in primes:
        compute_reduced_euler_characteristic_mod_k(p)