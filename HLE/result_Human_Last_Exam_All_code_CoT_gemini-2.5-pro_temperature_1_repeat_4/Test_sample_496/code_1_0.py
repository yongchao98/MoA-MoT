def solve_total_rank():
    """
    Calculates the total rank of the SO(4)-equivariant cohomology ring of
    SO(4) \ X for degrees up to 100.

    The Poincare series for the equivariant cohomology ring A is given by:
    P_A(t) = (1 + t^3) / (1 - t^4)^2
           = (1 + t^3) * sum_{j=0 to inf} (j+1)*t^(4j)
           = sum_{j=0 to inf} (j+1)*t^(4j) + sum_{j=0 to inf} (j+1)*t^(4j+3)

    The rank of A in degree k, denoted a_k, is the coefficient of t^k.
    - If k = 4j, a_k = j+1.
    - If k = 4j+3, a_k = j+1.
    - Otherwise, a_k = 0.

    We sum the ranks a_k for all k <= 100.
    """
    max_degree = 100

    # Part 1: Sum of ranks for degrees of the form 4j
    # The condition is 4j <= 100, so j <= 25.
    j_max1 = max_degree // 4
    sum1 = 0
    for j in range(j_max1 + 1):
        # The rank is j+1
        sum1 += (j + 1)

    # Part 2: Sum of ranks for degrees of the form 4j+3
    # The condition is 4j+3 <= 100, so 4j <= 97, which means j <= 24.
    j_max2 = (max_degree - 3) // 4
    sum2 = 0
    for j in range(j_max2 + 1):
        # The rank is j+1
        sum2 += (j + 1)

    total_rank = sum1 + sum2
    
    # Per instructions, outputting the numbers in the final equation
    print(f"Contribution from degrees of the form 4j (for j from 0 to {j_max1}): {sum1}")
    print(f"Contribution from degrees of the form 4j+3 (for j from 0 to {j_max2}): {sum2}")
    print(f"Final Equation: {sum1} + {sum2} = {total_rank}")
    print(f"The total rank of A in degree * <= 100 is: {total_rank}")

solve_total_rank()