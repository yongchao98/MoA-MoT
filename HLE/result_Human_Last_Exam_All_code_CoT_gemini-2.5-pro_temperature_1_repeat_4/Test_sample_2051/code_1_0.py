def solve_moduli_volume_degree():
    """
    This function calculates the degree of the polynomial Z_{g, n_+, n_-}
    for the specific case given in the problem.
    """
    # (a) Does the property of piecewise polynomiality of Z imply continuity?
    # Yes. The volume Z is an integral over a polytope in the space of edge lengths.
    # This polytope varies continuously with the boundary lengths L, so its volume
    # is also a continuous function of L. The polynomial pieces must match at
    # the boundaries between cells.
    answer_a = "Yes"

    # (b) Determine the degree of the polynomial Z_{0,3,1}.
    # The parameters for the specific case.
    g = 0
    n_plus = 3
    n_minus = 1

    # The degree depends on the genus g and the total number of boundaries n.
    n = n_plus + n_minus

    # The degree of the volume polynomial Z_{g,n} is given by the formula:
    # Degree = 6g - 6 + 2n.
    # This comes from the fact that the volume is calculated by integrating over
    # the E = 6g - 6 + 3n edge lengths of a trivalent fatgraph, subject to n
    # linear constraints from the boundary lengths. The degree is E - n.
    degree = 6 * g - 6 + 2 * n

    # Print the final formatted answer for both parts
    print(f"(a) {answer_a}")
    print(f"(b) To find the degree of Z_{{{g},{n_plus},{n_minus}}}, we first calculate the total number of boundaries, n:")
    print(f"    n = n_+ + n_- = {n_plus} + {n_minus} = {n}")
    print(f"    Then, we use the degree formula: Degree = 6*g - 6 + 2*n")
    print(f"    Degree = 6*({g}) - 6 + 2*({n}) = {6*g} - 6 + {2*n} = {degree}")

solve_moduli_volume_degree()
