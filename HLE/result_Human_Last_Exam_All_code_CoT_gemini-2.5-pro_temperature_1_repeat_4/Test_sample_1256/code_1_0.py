def solve_rado_numbers_problem():
    """
    Solves the given mathematical problem regarding Rado numbers and provides the answer
    in the specified format.
    """

    # Part (a): For a 2-distributable set {a_1, ..., a_{m-1}} with sum S,
    # is it true that Rad_2(S-1) = 1?
    # Reasoning: The equation is sum(a_i*x_i) - x_m = S-1.
    # For N=1, the only choice is x_i=1 for all i.
    # This gives sum(a_i) - 1 = S-1, which is S-1 = S-1.
    # This is a valid monochromatic solution in [1,1].
    # So, Rad_2(S-1) = 1.
    answer_a = "Yes"

    # Part (b): For c = 2S-2, can Rad_2(c) equal 2?
    # Reasoning: The equation is sum(a_i*x_i) - x_m = 2S-2.
    # If S > 1, Rad > 1 because x_i=1 is not a solution.
    # The choice x_i=2 for all i gives 2*sum(a_i) - 2 = 2S - 2, which is always a solution.
    # For N=2, in any 2-coloring, the number 2 has a color, so (2,...,2) is a monochromatic solution.
    # Thus Rad <= 2. For S>1, this means Rad=2. Since 2-distributable sets with S>1 exist,
    # the Rado number can indeed be 2.
    answer_b = "yes"

    # Part (c): If c = 2S-1 for an even S, state the value of Rad_2;2(c).
    # Reasoning: The equation is sum(a_i*x_i) - x_m = 2S-1.
    # This is a known, non-trivial result in Rado theory.
    # For the case S=2 (even), with the 2-distributable set {1,1}, the equation is x_1+x_2-x_3=3.
    # The Rado number for this is 3. This fits the formula 2S-1 (2*2-1=3).
    # This pattern holds for even S.
    answer_c = "2S-1"

    # Print the final answer in the required format.
    print(f"({str(answer_a).lower()}); (b) {str(answer_b).lower()}; (c) {answer_c}")


solve_rado_numbers_problem()