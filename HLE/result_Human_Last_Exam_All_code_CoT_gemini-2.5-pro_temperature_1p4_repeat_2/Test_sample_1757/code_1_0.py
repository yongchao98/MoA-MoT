import math

def solve_high_dimensional_ratio_sum():
    """
    This function solves the described problem by numerically summing the series of ratios.

    The problem asks for S = sum_{d=1 to infinity} (A_d / B_d), where:
    A_d: The expected Lebesgue measure (volume) of the convex hull of 
         {Origin, p_1, ..., p_d}, with p_i on the i-th axis.
    B_d: The expected pairwise Euclidean distance between any two points in that set.
    """

    # K is a constant representing the expected distance between two random points on
    # orthogonal axes within the unit ball: E[sqrt(t_1^2 + t_2^2)]
    # where t_1, t_2 are independent random variables from Uniform(-1, 1).
    # Its analytical value is (sqrt(2) + ln(1 + sqrt(2))) / 3.
    K = (math.sqrt(2) + math.log(1 + math.sqrt(2))) / 3.0

    total_sum = 0.0
    d = 1
    
    print("Calculating the sum of Ratio_d = A_d / B_d for d = 1, 2, 3, ...")
    print(f"Using the constant K = {K:.6f}\n")
    print("="*40)

    while True:
        # Calculate A_d = (1/d!) * (0.5)^d
        # This is the expected volume of the d-simplex.
        try:
            # We use log-gamma (lgamma) for better numerical stability with factorials.
            # lgamma(n+1) is log(n!).
            log_A_d = -math.lgamma(d + 1) - d * math.log(2)
            A_d = math.exp(log_A_d)
        except (ValueError, OverflowError):
            # This happens for very large d, but we will have converged by then.
            break

        # Calculate B_d = (1/(d+1)) * (1 + (d-1)*K)
        # This is the average of expected distances over all pairs in the set.
        B_d = (1.0 / (d + 1.0)) * (1.0 + (d - 1.0) * K)
        
        # Calculate the ratio for the current dimension.
        ratio_d = A_d / B_d
        
        # Output the numbers involved in the calculation for this dimension's term.
        print(f"Dimension d = {d}:")
        print(f"  A_d (Expected Volume)   = {A_d:.6e}")
        print(f"  B_d (Expected Distance)   = {B_d:.6f}")
        print(f"  Term to be added (Ratio)  = {ratio_d:.6f}")
        print("-" * 40)

        # Break the loop if the term becomes too small, indicating convergence.
        if ratio_d < 1e-15:
            break
            
        total_sum += ratio_d
        d += 1
    
    print("The summation has converged.")
    print(f"\nFinal Sum = {total_sum}")
    print(f"Result with three-decimal precision: {total_sum:.3f}")

solve_high_dimensional_ratio_sum()