import math

def solve():
    """
    Calculates the integer part of 10^4 * lim(F(N)/ln(N)) as N -> infinity.
    """
    # The problem reduces to finding two families of solutions corresponding to invariants k=13 and k=5.
    # The number of solutions F(N) for 1 <= a,b <= N is asymptotically given by:
    # F(N) ~ (2/ln(alpha_13) + 2/ln(alpha_5)) * ln(N)
    # where alpha_k is the largest root of x^2 - kx + 1 = 0.

    k1 = 13
    k2 = 5

    # Largest root for k=13 is (13 + sqrt(13^2 - 4))/2 = (13 + sqrt(165))/2
    val1_in_sqrt = k1**2 - 4
    alpha1 = (k1 + math.sqrt(val1_in_sqrt)) / 2

    # Largest root for k=5 is (5 + sqrt(5^2 - 4))/2 = (5 + sqrt(21))/2
    val2_in_sqrt = k2**2 - 4
    alpha2 = (k2 + math.sqrt(val2_in_sqrt)) / 2

    # The limit L = F(N)/ln(N) as N -> infinity
    limit_val = 2 / math.log(alpha1) + 2 / math.log(alpha2)

    # We need to find the integer part of 10^4 * L
    final_value = 10**4 * limit_val

    # As requested, printing the numbers in the final equation
    print(f"The expression for the limit is: 2/ln(({k1}+sqrt({val1_in_sqrt}))/2) + 2/ln(({k2}+sqrt({val2_in_sqrt}))/2)")
    print(f"The value to compute is: 10000 * ({limit_val})")
    print(f"The integer part of the result is: {int(final_value)}")

solve()