def solve_quiver_automorphism_problem():
    """
    This script determines the answers to the three questions based on a step-by-step
    deduction process as outlined in the thought plan.
    """

    # --- Step 1: Determine the answers for each question ---

    # Answer for (a): No.
    # The equivariance condition `g*sigma = lambda*sigma*g` leads to c_11*(lambda^2 - 1) = 0.
    # If lambda^2 = 1 (a common scenario, derived from (c) under unitarity assumptions),
    # this does not force the coefficient c_11 of a_j in sigma(a_j) to be zero.
    answer_a = "No"

    # Answer for (b): no.
    # The premises lead to the relation c_j^* = (lambda / mu_j) * mu_{j-1}^* * c_j.
    # For this to equal -mu_j^{-1} * c_j, it would require lambda * mu_{j-1}^* = -1.
    # This condition is not guaranteed to hold for arbitrary parameters lambda and mu.
    answer_b = "no"

    # Answer for (c): yes.
    # Assuming a typo in the question ("not intersected" should be "is intersected"),
    # the condition lambda^2 * mu_i * mu_i^* = 1 is necessary. This follows from the
    # g^2-equivariance of sigma (`g^2*sigma = lambda^2*sigma*g^2`) and the assumption that
    # g^2 preserves the algebra's relations (`g^2(sigma(a_i)) = sigma(a_i)`).
    answer_c = "yes"
    
    # --- Step 2: Print the final answer in the specified format ---

    # The problem asks for the output in a specific string format.
    print(f"(a) [{answer_a}]; (b) [{answer_b}]; (c) [{answer_c}].")

solve_quiver_automorphism_problem()