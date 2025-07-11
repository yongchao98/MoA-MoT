def solve_moduli_volume_question():
    """
    Solves the two-part question about the moduli space volume Z.
    """

    # Part (a): Continuity
    # The volume of the moduli space, Z, is a piecewise polynomial.
    # For a geometric quantity like volume, the different polynomial pieces
    # must meet continuously at the boundaries of their domains (cells).
    # A discontinuity would imply a non-physical jump in volume for an
    # infinitesimal change in boundary length.
    # Therefore, the function is continuous.
    answer_a = "Yes"

    # Part (b): Degree of the polynomial for g=0, n_+=3, n_-=1
    g = 0
    n_plus = 3
    n_minus = 1

    # The total number of boundaries is n.
    n = n_plus + n_minus

    # For the special case of genus g=0, the degree of the volume polynomial
    # Z_0,n is given by the formula: n - 2.
    if g == 0:
        degree_b = n - 2
    else:
        # For g > 0, the formula is 3g - 3 + n.
        degree_b = 3 * g - 3 + n

    # The question asks for the answer in a specific format.
    # We will print the final result to the console.
    final_answer_string = f"(a) [{answer_a}]; (b) [{degree_b}]."
    
    print("The final answer is formatted as: (a) [Yes/No]; (b) [Degree].")
    print(final_answer_string)
    
    # Also printing the required final answer format for the platform.
    print(f"<<<{final_answer_string}>>>")

solve_moduli_volume_question()