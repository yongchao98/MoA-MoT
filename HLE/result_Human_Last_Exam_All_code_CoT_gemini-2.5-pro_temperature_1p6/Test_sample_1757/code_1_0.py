import math

def solve_and_print_sum():
    """
    Calculates the sum over all natural dimensions of the specified ratio.

    The problem asks for S = sum_{d=1 to infinity} R_d, where R_d is a ratio.
    - Ratio: R_d = A_d / B_d
    - Numerator A_d: Expected volume of the convex hull of {O, p_1, ..., p_d}.
      - p_i = u_i * e_i, with u_i ~ U(-1, 1).
      - A_d = E[ (1/d!) * |det(p_1, ..., p_d)| ] = (1/d!) * (E[|u|])^d = (1/d!) * (1/2)^d.
    - Denominator B_d: Expected pairwise distance between points in {O, p_1, ..., p_d}.
      - This is the average of C(d+1, 2) pairwise distance expectations.
      - E[dist(O, p_i)] = E[|u_i|] = 1/2.
      - E[dist(p_i, p_j)] = E[sqrt(u_i^2 + u_j^2)]. Let's call this C_dist.
        C_dist = (1/3) * (sqrt(2) + log(1 + sqrt(2))).
      - After simplification, B_d = (1 + (d-1) * C_dist) / (d+1).
    - Final Ratio R_d = A_d / B_d = (d+1) / (d! * 2^d * (1 + (d-1) * C_dist)).
    """
    # Define the constant for the expected distance between two points on orthogonal axes.
    # C_dist = E[sqrt(x^2 + y^2)] for x, y ~ U(-1, 1).
    try:
        C_dist = (math.sqrt(2) + math.log(1 + math.sqrt(2))) / 3
    except ImportError:
        # Fallback for systems without math.log etc.
        C_dist = 0.7651957

    total_sum = 0.0
    
    print("This script calculates the sum S = R_1 + R_2 + R_3 + ...")
    print("The term R_d for each dimension d is the ratio of an expected volume to an expected distance.")
    print("The formula derived for the term is: R_d = (d+1) / (d! * 2^d * (1 + (d-1) * C_dist))\n")
    print("--- Calculating the terms of the sum ---")

    # Loop through dimensions d, summing the terms R_d until they are negligible.
    for d in range(1, 30): 
        try:
            # Numerator of the final expression for R_d
            numerator_r = d + 1
            
            # Components of the denominator
            factorial_d = math.factorial(d)
            power_of_2_d = 2**d
            dist_term = 1 + (d - 1) * C_dist
            
            # Denominator of the final expression for R_d
            denominator_r = factorial_d * power_of_2_d * dist_term
            
            # The term R_d for the current dimension
            R_d = numerator_r / denominator_r

            # "output each number in the final equation"
            # Here we print the breakdown for each R_d which constitutes the sum.
            print(f"Term for d = {d}:")
            print(f"  R_{d} = ({d}+1) / ({d}! * 2^{d} * (1 + ({d}-1) * {C_dist:.4f}))")
            print(f"      = {numerator_r} / ({factorial_d} * {power_of_2_d} * {dist_term:.4f})")
            print(f"      = {numerator_r} / {denominator_r:.4f}")
            print(f"      ≈ {R_d:.8f}\n")

            total_sum += R_d

            # If a term is very small, subsequent terms will be even smaller.
            # Stop summing as they won't affect the result at the required precision.
            if R_d < 1e-12:
                print(f"Term R_{d} is smaller than 1e-12, stopping the summation.")
                break
        except OverflowError:
            print(f"OverflowError at d={d}. The sum has effectively converged.")
            break
    
    print("--- Final Result ---")
    print(f"The total sum S converges to: {total_sum}")
    print(f"The sum rounded to three-decimal precision is: {total_sum:.3f}")

if __name__ == '__main__':
    solve_and_print_sum()
