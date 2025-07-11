def solve_smallest_n():
    """
    This function finds the smallest integer n >= 2 based on the logic derived from the problem statement.
    It calculates the smallest possible n for each valid case and then determines the overall minimum.
    """

    print("Finding the smallest integer n >= 2 with the specified properties.")
    print("This involves analyzing three possible cases for n based on modular arithmetic.\n")

    # Case B: n is even and n = 1 (mod 5^9)
    # n = k * 5^9 + 1. n must be even, so k must be odd. Smallest positive k is 1.
    power_of_5 = 5**9
    n_B = 1 * power_of_5 + 1
    print("Candidate from Case B (n = 0 mod 2, n = 1 mod 5^9):")
    print(f"n = 1 * 5^9 + 1 = 1 * {power_of_5} + 1 = {n_B}\n")

    # Case C: n = 1 (mod 2^9) and n is a multiple of 5
    # This also requires n != 1 (mod 2^10).
    # n = k * 2^9 + 1 = 512*k + 1. n must be a multiple of 5.
    # 512*k + 1 = 0 (mod 5) => 2*k + 1 = 0 (mod 5) => k = 2 (mod 5).
    # k=2 gives n=1025, which is 1 (mod 1024), so it's invalid.
    # The next k is 7.
    power_of_2 = 2**9
    n_C = 7 * power_of_2 + 1
    print("Candidate from Case C (n = 1 mod 2^9, n = 0 mod 5, n != 1 mod 2^10):")
    print(f"n = 7 * 2^9 + 1 = 7 * {power_of_2} + 1 = {n_C}\n")

    # Case D: n = 1 (mod 10^9)
    # This also requires n != 1 (mod 10^10).
    # n = k * 10^9 + 1. k cannot be a multiple of 10. Smallest positive k is 1.
    power_of_10 = 10**9
    n_D = 1 * power_of_10 + 1
    print("Candidate from Case D (n = 1 mod 10^9, n != 1 mod 10^10):")
    print(f"n = 1 * 10^9 + 1 = 1 * {power_of_10} + 1 = {n_D}\n")

    # The result is the minimum of the candidates from the valid cases.
    result = min(n_B, n_C, n_D)

    print("The smallest integer is the minimum of these candidates.")
    print(f"Final Answer: min({n_B}, {n_C}, {n_D}) = {result}")

solve_smallest_n()