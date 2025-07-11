def solve_moduli_volume_questions():
    """
    Solves the two questions about the volume of the moduli space Z.
    """
    
    # Part (a): On the continuity of Z.
    # The volume of a moduli space, Z, is a geometric quantity. While it is
    # piecewise polynomial, the polynomial pieces join continuously at the
    # boundaries of the cells. Thus, the function is continuous.
    answer_a = "Yes"

    # Part (b): Determine the degree of Z_{0,3,1}.
    # The parameters are given as:
    g = 0
    n_plus = 3
    n_minus = 1

    # The total number of boundaries is n = n_+ + n_-.
    n_total = n_plus + n_minus
    
    # The degree of the polynomial Z_{g, n_+, n_-} is given by the formula:
    # degree = 3g - 3 + n
    degree = 3 * g - 3 + n_total

    # As requested, we will output the equation used for the calculation.
    print(f"Degree calculation: 3*g - 3 + (n_+ + n_-) = 3*{g} - 3 + ({n_plus} + {n_minus}) = {degree}")

    # Print the final answer in the required format.
    print(f"(a) {answer_a}; (b) {degree}")

solve_moduli_volume_questions()