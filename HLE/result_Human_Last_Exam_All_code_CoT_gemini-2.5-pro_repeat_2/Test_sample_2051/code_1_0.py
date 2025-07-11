def solve_moduli_volume_problem():
    """
    Solves the two-part question about the volume of the moduli space of oriented metric ribbon graphs.
    """

    # Part (a): Continuity
    # The property of being "piecewise polynomial" does not, in general, imply continuity.
    # However, the function Z_{g, n_+, n_-} represents the volume of a moduli space.
    # As a geometric quantity, this volume is known to be a continuous function of the
    # boundary length parameters. The different polynomial pieces join together seamlessly
    # at the boundaries between cells.
    answer_a = "Yes"

    # Part (b): Degree of the polynomial
    # We are asked to find the degree of Z_{g, n_+, n_-} for g = 0, n_+ = 3, and n_- = 1.
    g = 0
    n_plus = 3
    n_minus = 1

    # From the theory of intersection numbers on moduli spaces and the work of P. Norbury,
    # the degree of the polynomial Z_{g, n_+, n_-} is given by the formula:
    # degree = 2g - 2 + n_+ + 2n_-
    # We will use this formula to calculate the degree.

    print("Calculation for the degree of the polynomial Z_{0,3,1}:")
    print("---------------------------------------------------------")
    print("The formula for the degree is: d = 2g - 2 + n_+ + 2n_-")
    print(f"Given parameters are: g = {g}, n_+ = {n_plus}, n_- = {n_minus}")
    
    print("\nStep 1: Substitute the parameters into the formula.")
    print(f"d = 2*({g}) - 2 + {n_plus} + 2*({n_minus})")

    # Perform the calculation step-by-step
    term1 = 2 * g
    term4 = 2 * n_minus
    degree = term1 - 2 + n_plus + term4
    
    print("\nStep 2: Evaluate the terms and simplify.")
    print(f"d = {term1} - 2 + {n_plus} + {term4}")
    print(f"d = {term1 - 2} + {n_plus} + {term4}")
    print(f"d = {term1 - 2 + n_plus} + {term4}")
    print(f"d = {degree}")
    
    answer_b = degree
    
    print("\n---------------------------------------------------------")
    print("Final Answer:")
    print(f"(a) {answer_a}; (b) {answer_b}")

solve_moduli_volume_problem()