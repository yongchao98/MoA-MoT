def check_coadjoint_orbit_properties():
    """
    Provides answers and explanations to the questions on coadjoint orbits.
    """
    
    # --- Question (a) ---
    answer_a = "True"

    # --- Question (b) ---
    answer_b = "No"
    # Set up the counterexample for SU(n)
    n = 3
    # The formula from the question
    betti_formula_result = n - 1
    # The actual Betti number for a singular orbit
    # The orbit O_lambda for a singular lambda can be CP^(n-1), e.g., CP^2 for SU(3)
    actual_betti_number = 1
    
    # --- Question (c) ---
    answer_c = "No"

    print("Answers to the questions on coadjoint orbits:\n")
    
    print(f"(a) {answer_a}. Every coadjoint orbit of a compact semisimple Lie group is a homogeneous Kahler manifold.")
    
    print(f"\n(b) {answer_b}. The statement is only true for regular lambda. For singular lambda, it fails.")
    print("    Counterexample: For G = SU(n), let n = {0}.".format(n))
    print("    The rank is n - 1, so the value from the question's formula is {0} - 1 = {1}.".format(n, betti_formula_result))
    print("    However, for a singular lambda on a wall of the Weyl alcove, the orbit can be the complex projective space CP^{0},".format(n-1))
    print("    which has a second Betti number b2 = {0}.".format(actual_betti_number))
    print("    Final comparison: The actual b2({0}) is not equal to the predicted {1}.".format(actual_betti_number, betti_formula_result))
    print("    The final equation is: {0} != {1}".format(actual_betti_number, betti_formula_result))

    print(f"\n(c) {answer_c}. The GKM conditions (e.g., pairwise linear independence of tangential weights at fixed points) are not satisfied for all coadjoint orbits, with known counterexamples for partial flag manifolds in non-simply-laced Lie groups.")

check_coadjoint_orbit_properties()