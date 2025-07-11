import math

def solve_pm(m):
    """
    Calculates the probability P_m for a given positive integer m.

    The problem is to find the probability that a sequence a_1, ..., a_{4m+2}
    is (i,j)-divisible. This has been shown to be equivalent to finding the
    number of pairs (i, j) that can be removed from {1, ..., 4m+2} such that
    the remaining 4m indices can be partitioned into m arithmetic progressions of length 4.

    The number of favorable pairs (i,j) has been determined to be m*(m+5)/2.
    The total number of pairs (i,j) is C(4m+2, 2).
    The probability P_m is the ratio of these two quantities.
    """
    if not isinstance(m, int) or m <= 0:
        print("Error: m must be a positive integer.")
        return

    # Numerator of the probability P_m
    # Number of favorable pairs = m*(m+5)/2
    # We will handle the division by 2 at the end
    favorable_numerator = m * (m + 5)
    favorable_denominator = 2

    # Denominator of the probability P_m
    # Total number of pairs = C(4m+2, 2) = (4m+2)*(4m+1)/2
    total_numerator = (4 * m + 2) * (4 * m + 1)
    total_denominator = 2
    
    # P_m = (favorable_numerator / favorable_denominator) / (total_numerator / total_denominator)
    # P_m = favorable_numerator / total_numerator
    
    final_numerator = favorable_numerator
    final_denominator = total_numerator
    
    # Simplify the fraction by dividing by the greatest common divisor
    common_divisor = math.gcd(final_numerator, final_denominator)
    
    simplified_numerator = final_numerator // common_divisor
    simplified_denominator = final_denominator // common_divisor
    
    print(f"For m = {m}:")
    print(f"The number of favorable pairs (i, j) is: {favorable_numerator // favorable_denominator}")
    print(f"The total number of pairs (i, j) is: {total_numerator // total_denominator}")
    
    # Final equation format requested in the problem.
    # We will show the components of the derived formula for Pm
    # P_m = m(m+5) / (2 * (2m+1) * (4m+1))
    
    term1_num = m
    term2_num = m + 5
    
    term1_den = 2
    term2_den = 2 * m + 1
    term3_den = 4 * m + 1
    
    calc_num = term1_num * term2_num
    calc_den = term1_den * term2_den * term3_den
    
    common = math.gcd(calc_num, calc_den)

    print("\nThe probability P_m can be calculated with the formula: m * (m+5) / (2 * (2*m+1) * (4*m+1))")
    print(f"Final Equation: P_m = ({term1_num} * {term2_num}) / ({term1_den} * {term2_den} * {term3_den})")
    print(f"Calculated Probability, P_{m}: {calc_num//common} / {calc_den//common}")


# Example usage for m=1
solve_pm(1)

# Example usage for m=2
solve_pm(2)