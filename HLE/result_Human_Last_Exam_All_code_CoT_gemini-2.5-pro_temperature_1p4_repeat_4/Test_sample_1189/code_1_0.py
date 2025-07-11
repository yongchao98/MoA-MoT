import math
import sys

def combinations(n, k):
    """
    Calculates the number of combinations (n choose k).
    Returns 0 if k > n or k < 0.
    """
    if k < 0 or k > n:
        return 0
    if k == 0 or k == n:
        return 1
    if k > n // 2:
        k = n - k
    
    res = 1
    for i in range(k):
        res = res * (n - i) // (i + 1)
    return res

def solve():
    """
    Solves the problem for given n and m.
    Reads n and m from command line arguments.
    """
    if len(sys.argv) != 3:
        print("Usage: python your_script_name.py <n> <m>")
        return

    try:
        n = int(sys.argv[1])
        m = int(sys.argv[2])
        if n <= 0 or m <= 0:
            raise ValueError("n and m must be positive integers.")
    except ValueError as e:
        print(f"Error: {e}")
        return

    # Term 1: C(2^n - 1, m)
    t1_n = 2**n - 1
    t1_k = m
    t1_val = combinations(t1_n, t1_k)

    # Term 2: (2^n - 1)
    t2_val = 2**n - 1

    # Term 3: (-1)^ceil(m/2)
    t3_val = 1 if math.ceil(m / 2) % 2 == 0 else -1

    # Term 4: C(2^(n-1) - 1, floor(m/2))
    t4_n = 2**(n - 1) - 1
    t4_k = m // 2
    t4_val = combinations(t4_n, t4_k)
    
    # Calculate final answer
    numerator = t1_val + t2_val * t3_val * t4_val
    denominator = 2**n
    final_answer = numerator // denominator

    # Print the equation breakdown
    print(f"For n = {n} and m = {m}, the formula is:")
    print(f"f({n}, {m}) = ( C({t1_n}, {t1_k}) + ({t2_val}) * ({t3_val}) * C({t4_n}, {t4_k}) ) / {denominator}")
    print(f"       = ( {t1_val} + {t2_val} * {t3_val} * {t4_val} ) / {denominator}")
    print(f"       = ( {t1_val} + {t2_val * t3_val * t4_val} ) / {denominator}")
    print(f"       = {numerator} / {denominator}")
    print(f"       = {final_answer}")
    
    print(f"\n<<< {final_answer} >>>")


if __name__ == '__main__':
    solve()
