import math

def solve():
    """
    Calculates the number of 1324-avoiding permutations of length n with k inversions,
    for n=333 and k=3.
    """
    n = 333
    k = 3

    # The problem is to find av_n^k(1324) for n=333, k=3.
    # A known formula for n >= 5 is:
    # av_n^3(1324) = C(n-2, 3) + C(n-2, 2) + C(n-2, 1) + 2
    
    if n < 5:
        print(f"The formula is valid for n >= 5. n={n} is not supported by this formula.")
        return

    m = n - 2
    
    # Calculate each term in the formula
    try:
        # math.comb(n, k) calculates the binomial coefficient "n choose k"
        c3 = math.comb(m, 3)
        c2 = math.comb(m, 2)
        c1 = math.comb(m, 1)
        constant = 2
    except ValueError as e:
        print(f"Error calculating combinations: {e}")
        return

    # Sum the terms to get the final result
    result = c3 + c2 + c1 + constant

    print(f"To find the value of av_{n}^{k}(1324), we use the formula:")
    print(f"C(n-2, 3) + C(n-2, 2) + C(n-2, 1) + 2")
    print(f"\nFor n = {n}, we have n-2 = {m}. The equation becomes:")
    print(f"C({m}, 3) + C({m}, 2) + C({m}, 1) + 2")
    print(f"\nCalculating each term:")
    print(f"C({m}, 3) = {c3}")
    print(f"C({m}, 2) = {c2}")
    print(f"C({m}, 1) = {c1}")
    print(f"\nSumming the terms:")
    print(f"{c3} + {c2} + {c1} + {constant} = {result}")

solve()