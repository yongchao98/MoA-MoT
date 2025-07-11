import math

def solve():
    """
    Calculates the sum described in the problem statement.
    """
    # Constant C_sqrt = E[sqrt(t_i^2 + t_j^2)]
    C_sqrt = (math.sqrt(2) + math.log(1 + math.sqrt(2))) / 3.0

    total_sum = 0.0
    
    # The equation will be built by summing terms R_d
    print("The sum is calculated as S = R(1) + R(2) + R(3) + ...")
    print("The terms of the sum are:")

    # Loop over dimensions d. The terms decrease very rapidly, so a small range is sufficient.
    for d in range(1, 10):
        # A(d) = 1 / (2d)^d
        try:
            A_d = 1.0 / math.pow(2 * d, d)
        except OverflowError:
            # For large d, A_d becomes effectively zero
            A_d = 0.0

        # B(d) is the expected pairwise distance
        if d == 1:
            B_d = 0.5
        else:
            # B_d = (5d - 2 + 3(d-1)^2 * C_sqrt) / (3d(d+1))
            numerator = 5 * d - 2 + 3 * math.pow(d - 1, 2) * C_sqrt
            denominator = 3 * d * (d + 1)
            B_d = numerator / denominator

        if B_d == 0:
            R_d = 0.0
        else:
            R_d = A_d / B_d
        
        # Output each number in the final equation
        print(f"R({d}) = {R_d:.7f}")
        
        total_sum += R_d

        # Break if the term is too small to affect the result to 3 decimal places
        if R_d < 1e-5:
            break
            
    print("\nFinal Result:")
    # Print the sum formatted to three decimal places
    print(f"The sum S is approximately: {total_sum:.3f}")

solve()
<<<1.117>>>