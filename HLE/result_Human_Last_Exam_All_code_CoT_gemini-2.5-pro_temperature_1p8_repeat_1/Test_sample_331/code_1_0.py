def solve_chow_group_problem():
    """
    This function stores and prints the computed pairs (m(X), M(X)) for the four given varieties.
    The values are derived from known results in algebraic geometry.
    """
    
    # For X_1, a genus 2 curve: m(X_1) = 2, M(X_1) = 3
    # m=2 is achieved at the 6 Weierstrass points.
    # M=3 is the value for any other point.
    x1_m = 2
    x1_M = 3
    pair1 = (x1_m, x1_M)
    
    # For X_2, a general genus 7 curve: m(X_2) = 7, M(X_2) = 7
    # For a general curve of genus g, edeg(C, p) = g for all points p.
    x2_m = 7
    x2_M = 7
    pair2 = (x2_m, x2_M)

    # For X_3, an Enriques surface: m(X_3) = 1, M(X_3) = 1
    # For an Enriques surface S, CH_0(S) is isomorphic to Z, so [p]=[q] for all p,q.
    x3_m = 1
    x3_M = 1
    pair3 = (x3_m, x3_M)
    
    # For X_4, the Grassmannian G(3,6): m(X_4) = 1, M(X_4) = 1
    # G(3,6) is a rationally connected variety, so CH_0(G(3,6)) is isomorphic to Z.
    x4_m = 1
    x4_M = 1
    pair4 = (x4_m, x4_M)

    # Formatting the output as requested.
    # The output string should be like: (2, 3), (7, 7), (1, 1), (1, 1)
    # The parentheses and commas are part of the formatting.
    # Note that print(pair) automatically adds the parentheses.
    print(f"({pair1[0]}, {pair1[1]}), ({pair2[0]}, {pair2[1]}), ({pair3[0]}, {pair3[1]}), ({pair4[0]}, {pair4[1]})")

solve_chow_group_problem()