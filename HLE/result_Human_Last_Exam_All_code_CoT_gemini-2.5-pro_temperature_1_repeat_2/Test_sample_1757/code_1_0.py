import math

def solve_sum():
    """
    Calculates the sum over all natural dimensions of the specified ratio.

    The plan is as follows:
    1. Calculate the constant K, which is the expected distance between two points
       on orthogonal axes within the unit cube, i.e., E[sqrt(t1^2 + t2^2)]
       for t1, t2 ~ U(-1, 1).
    2. Initialize a total sum to 0.0.
    3. Loop through dimensions d = 1, 2, 3, ...
    4. In each iteration, calculate the term R_d = A_d / B_d, where:
       - A_d = (1/d!) * (1/2)^d is the expected volume.
       - B_d = (1 + (d-1)*K) / (d+1) is the expected pairwise distance.
    5. Add R_d to the total sum.
    6. The loop continues until the term R_d is smaller than a precision
       threshold (e.g., 1e-15), ensuring the sum converges.
    7. Print the final sum formatted to three decimal places.
    """
    # Step 1: Calculate the constant K
    # K = (sqrt(2) + ln(1+sqrt(2))) / 3. Note that ln(1+sqrt(2)) = asinh(1).
    K = (math.sqrt(2) + math.asinh(1)) / 3.0

    # Step 2: Initialize sum
    total_sum = 0.0

    # Step 3-6: Loop through dimensions and sum the terms R_d
    # The loop will be stopped when the term added is too small to affect the result.
    # A max dimension of 100 is more than sufficient.
    for d in range(1, 101):
        try:
            # Numerator of R_d
            numerator_Rd = float(d + 1)

            # Denominator of R_d
            # Using math.pow for 2**d to ensure float arithmetic
            denominator_Rd = math.factorial(d) * math.pow(2, d) * (1 + (d - 1) * K)

            # Calculate R_d
            r_d = numerator_Rd / denominator_Rd
            
            # Add to sum
            total_sum += r_d

            # Convergence check
            if r_d < 1e-15:
                break
        except OverflowError:
            # math.factorial(d) will grow very large.
            # If an overflow occurs, the term is certainly negligible.
            break

    # Step 7: Print the final result formatted as requested.
    # The equation is S = [Sum of R_d], we print the number S.
    print(f"{total_sum:.3f}")

solve_sum()