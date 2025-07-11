def solve_moduli_volume_questions():
    """
    This function addresses the user's questions about the properties
    of the moduli space volume Z_{g, n_+, n_-}.
    """

    # Part (a): Continuity of the piecewise polynomial function Z.
    # The volume of a moduli space is a continuous function of its parameters.
    # Although Z is piecewise polynomial, the polynomial pieces join continuously
    # at the boundaries of the cells in the parameter space.
    answer_a = "Yes"

    # Part (b): Degree of the polynomial Z_{0,3,1}.
    # The parameters are g=0, n_+=3, n_-=1.
    g = 0
    n_plus = 3
    n_minus = 1

    # The total number of boundaries n is n_+ + n_-.
    n = n_plus + n_minus

    # The degree of the volume polynomial Z_{g,n} is given by the
    # dimension of the moduli space of curves, which is 3g - 3 + n.
    degree = 3 * g - 3 + n
    answer_b = degree
    
    # Print the answers. The second part includes the calculation as requested.
    print(f"(a) {answer_a}")
    print(f"(b) The degree of the polynomial is determined by the formula 3g - 3 + n.")
    print(f"For g = {g}, n = {n}, the calculation is: 3*({g}) - 3 + {n} = {degree}")

solve_moduli_volume_questions()