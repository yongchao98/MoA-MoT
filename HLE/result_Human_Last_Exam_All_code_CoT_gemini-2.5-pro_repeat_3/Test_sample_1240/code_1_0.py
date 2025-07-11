def solve_lattice_questions():
    """
    Solves three questions about root systems of d-neighbors of Z^n.

    The method is based on the analysis of the condition u.v = 0 (mod d)
    for a root v to be in the visible part R_2(M).
    """

    # --- Part (a) ---
    # Question: For a d-neighbor N of Z^12, can R_2(M) be of type A_11?
    # Plan: We need to find if there exist an integer d > 1 and a primitive vector
    # u in Z^12 such that R_2(M) = {v in D_12 | u.v = 0 (mod d)} is an A_11 system.
    # The roots of A_11 are {+/- (e_i - e_j) | 1 <= i < j <= 12}.
    # This requires:
    # 1. u . (e_i - e_j) = 0 (mod d) for all i, j.
    # 2. u . (e_i + e_j) != 0 (mod d) for all i, j.

    # Let's test u = (1, 1, ..., 1). This vector is primitive.
    # Condition 1: u . (e_i - e_j) = u_i - u_j = 1 - 1 = 0.
    # This is congruent to 0 (mod d) for any d. So this condition is met.
    # Condition 2: u . (e_i + e_j) = u_i + u_j = 1 + 1 = 2.
    # We need 2 != 0 (mod d), which means d must not be a divisor of 2.
    # We can choose d = 3. This satisfies the condition.
    # With u = (1, ..., 1) and d = 3, R_2(M) is precisely the A_11 root system.
    answer_a = "Yes"

    # --- Part (b) ---
    # Question: Can the visible root system R_2(M) of a d-neighbor N of Z^15
    # contain a D_7 component?
    # Plan: We need to find if there exists (d, u) such that R_2(M) contains a D_7
    # component. Let the component be on the index set S = {1, ..., 7}.
    # This requires:
    # 1. For i,j in S, +/-e_i +/- e_j must be in M. This means u_i +/- u_j = 0 (mod d).
    #    This implies u_i = u_j (mod d) and 2*u_i = 0 (mod d) for all i,j in S.
    # 2. For the component to be separate, for i in S and k not in S,
    #    +/-e_i +/- e_k must NOT be in M. This means u_i +/- u_k != 0 (mod d).

    # Let's try to construct such a (d, u). Let d = 3.
    # From 2*u_i = 0 (mod 3), since gcd(2, 3) = 1, it implies u_i = 0 (mod 3) for i in S.
    # For condition 2, for k not in S, we need u_i +/- u_k = 0 +/- u_k != 0 (mod 3).
    # This means u_k must not be a multiple of 3. Let's set u_k = 1 (mod 3).
    # We can choose a primitive vector u satisfying these conditions, for example:
    # u = (3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1).
    # This u is primitive because gcd(3, 1) = 1.
    # This construction works, and R_2(M) contains D_7 as a component.
    answer_b = "yes"

    # --- Part (c) ---
    # Question: For n = 18 and d = 5, is it possible for R_2(M) to include
    # more than one D_k component?
    # Plan: Assume R_2(M) contains two disjoint components, D_k1 on index set S1
    # and D_k2 on S2. We will check if this leads to a contradiction.
    # For D_k1 on S1, we need u_i = u_j (mod 5) and 2*u_i = 0 (mod 5) for i, j in S1.
    # Since gcd(2, 5) = 1, the second congruence 2*u_i = 0 becomes u_i = 0 (mod 5).
    # So, u_i = 0 (mod 5) for all i in S1.
    # Similarly, for D_k2 on S2, we must have u_j = 0 (mod 5) for all j in S2.
    # Now, for the components to be separate, roots connecting S1 and S2 must not be in M.
    # Let's check a connecting root v = e_i + e_j where i is in S1, j is in S2.
    # u . v = u_i + u_j.
    # Since u_i = 0 (mod 5) and u_j = 0 (mod 5), we have u_i + u_j = 0 + 0 = 0 (mod 5).
    # This means the connecting root v IS in M.
    # This contradicts the assumption that S1 and S2 form separate components.
    # Any two such D-type systems would merge into a single, larger D-type system.
    # Therefore, it is not possible.
    answer_c = "no"

    print(f"(a) [{answer_a}]; (b) [{answer_b}]; (c) [{answer_c}].")

solve_lattice_questions()
<<< (a) [Yes]; (b) [yes]; (c) [no]. >>>