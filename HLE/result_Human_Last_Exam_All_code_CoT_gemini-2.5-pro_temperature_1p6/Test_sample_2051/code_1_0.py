def solve_moduli_volume_questions():
    """
    This script addresses two questions about the volume of moduli spaces of ribbon graphs.
    (a) On the implication of continuity from piecewise polynomiality.
    (b) On the degree of a specific volume polynomial Z_{0,3,1}.
    """

    # --- Part (a): Continuity ---
    # The property of a function being "piecewise polynomial" means its domain can be
    # partitioned into a finite number of regions, and on each region, the function
    # is given by a polynomial. This property alone does not imply continuity.
    # A standard counterexample is a step function, which is composed of constant
    # (degree-0) polynomial pieces but is discontinuous at the points where it jumps.
    # Although the volume function Z is in fact continuous, this is a deeper property
    # that comes from the geometry of the moduli space, not from the piecewise
    # polynomial description alone.
    answer_a = "No"

    # --- Part (b): Degree Calculation ---
    # We are given the case g = 0, n_+ = 3, n_- = 1.
    g = 0
    n_plus = 3
    n_minus = 1

    # The degree of the polynomial Z_{g, n} is related to the (complex) dimension
    # of the corresponding moduli space M_{g,n}.

    # Step 1: Calculate the total number of boundaries, n.
    n = n_plus + n_minus

    # Step 2: Use the formula for the degree of the volume polynomial.
    # The volume Z_{g,n} is a polynomial in the squares of the boundary lengths (L_i^2).
    # The degree of Z as a polynomial in the variables {L_1^2, ..., L_n^2} is given
    # by the formula: d = 3g - 3 + n.
    degree_in_L_squared = 3 * g - 3 + n

    # Step 3: Convert the degree in L_i^2 to the degree in L_i.
    # If a function is a polynomial of degree 'd' in x^2, it is a polynomial of
    # degree 2d in x.
    final_degree = 2 * degree_in_L_squared
    answer_b = final_degree

    # --- Print the results ---
    print(f"(a) {answer_a}")
    print("-" * 20)
    print(f"(b) {answer_b}")
    print("Calculation Steps:")
    print(f"The parameters are g = {g}, n_+ = {n_plus}, n_- = {n_minus}.")
    print(f"The total number of boundaries is n = n_+ + n_- = {n_plus} + {n_minus} = {n}.")
    print(f"The volume is a polynomial in L_i^2 of degree d = 3g - 3 + n.")
    print(f"  d = 3*({g}) - 3 + {n} = {0 - 3 + n} = {degree_in_L_squared}")
    print(f"The degree of the polynomial in the lengths L_i is 2 * d.")
    print(f"  Degree = 2 * {degree_in_L_squared} = {final_degree}")

solve_moduli_volume_questions()