def solve_problem():
    """
    This function solves the problem by determining the cardinalities of S/A, B/S, and H_1(B/A, Q).
    """

    # Part 1: Determine the cardinality of S/A
    # The initial object A is the scale id_Z: Z -> Z. So, the group G_A is Z.
    # The scale S is the inclusion map S: Z -> *R (the hyperreals).
    # The canonical map from A to S is S itself, so its image is the standard integers Z_st within *R.
    # The quotient S/A is *R / Z_st.
    # The cardinality of the hyperreals, |*R|, is c = 2^aleph_0 = Beth_1.
    # The cardinality of the integers, |Z_st|, is aleph_0 = Beth_0.
    # From |*R| = |*R / Z_st| * |Z_st|, we have Beth_1 = |*R / Z_st| * Beth_0.
    # This implies that the cardinality of the quotient is Beth_1.
    cardinality_S_div_A = "Beth_1"

    # Part 2: Determine the cardinality of B/S
    # The terminal object B is the zero scale z_0: Z -> {0}. The group G_B is the trivial group {0}.
    # The canonical map from S (*R) to B ({0}) is the zero homomorphism.
    # The image of this map is the subgroup {0} within G_B.
    # The quotient B/S is G_B / im(h_SB) = {0} / {0}, which is the trivial group.
    # The cardinality of the trivial group is 1.
    cardinality_B_div_S = 1

    # Part 3: Determine the cardinality of H_1(B/A, Q)
    # The space B/A is the quotient G_B / im(h_AB).
    # The canonical map from A (Z) to B ({0}) is the zero homomorphism.
    # Its image is {0}.
    # The space B/A is {0} / {0}, which is the trivial group. As a topological space, it's a single point.
    # We need to find the cardinality of H_1({pt}, Q).
    # The homology of a point with integer coefficients is H_n({pt}, Z) = {0} for n > 0.
    # By the Universal Coefficient Theorem, H_1({pt}, Q) is also the trivial group {0}.
    # The cardinality of the trivial group is 1.
    cardinality_H1_B_div_A = 1

    # Print the final answer in the required format
    print(f"{cardinality_S_div_A} {cardinality_B_div_S} {cardinality_H1_B_div_A}")

solve_problem()