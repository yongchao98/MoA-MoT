import numpy as np

def solve_ququint_problem():
    """
    Calculates the most probable final state of the entangled ququint system
    after applying gate Q and measurement.
    """
    # 1. Define the gate Q as a matrix.
    # The problem defines Q|j> which corresponds to the j-th column of the matrix.
    # Q|0> = 1/sqrt(2) * (|1> + |2>) -> col 0
    # Q|1> = 1/sqrt(2) * (|0> + |3>) -> col 1
    # Q|2> = 1/sqrt(2) * (|1> + |4>) -> col 2
    # Q|3> = 1/sqrt(2) * (|2> + |0>) -> col 3
    # Q|4> = 1/sqrt(2) * (|3> + |2>) -> col 4
    Q = (1 / np.sqrt(2)) * np.array([
        [0, 1, 0, 1, 0],
        [1, 0, 1, 0, 0],
        [1, 0, 0, 1, 1],
        [0, 1, 0, 0, 1],
        [0, 0, 1, 0, 0]
    ])

    # 2. Determine the state of ququint B for each measurement outcome of A.
    # If we measure ququint A in state |j>, the unnormalized state of B is
    # given by the coefficients in the j-th row of the Q matrix.
    
    # 3. Calculate the probabilities for each measurement outcome.
    # The probability P(A=j) is proportional to the squared norm of the j-th row vector.
    # The total probability is normalized by Tr(Q_dagger @ Q) / d, which is 5/5 = 1.
    # So, P(j) = (1/5) * ||row_j(Q)||^2 * 5 = ||row_j(Q)||^2
    # Let's use the matrix Q directly. The state of B is proportional to the j-th row.
    
    probs = np.zeros(5)
    for j in range(5):
        # The vector for the state of B is the j-th row of Q
        b_state_vector = Q[j, :]
        # The probability is the squared norm of this vector
        probs[j] = np.linalg.norm(b_state_vector)**2

    # 4. Find the most probable outcome.
    most_probable_index = np.argmax(probs)
    
    # 5. Construct the final state for this outcome.
    # The state of A is the basis state corresponding to the most probable index.
    a_state_component = f"|{most_probable_index}>_A"

    # The state of B is given by the corresponding row vector, which needs to be normalized.
    # We use the raw integer coefficients for clarity before normalization.
    b_unnormalized_coeffs = Q[most_probable_index, :] * np.sqrt(2)
    
    # Calculate the normalization factor from the integer coefficients.
    norm_squared = int(round(np.sum(b_unnormalized_coeffs**2)))

    # Build the string for the state of B.
    b_state_terms = []
    for i, coeff in enumerate(b_unnormalized_coeffs):
        if not np.isclose(coeff, 0):
            b_state_terms.append(f"|{i}>_B")
    
    b_state_superposition = " + ".join(b_state_terms)
    b_state_component = f"(1/sqrt({norm_squared})) * ({b_state_superposition})"

    # 6. Print the final equation.
    final_equation = f"{a_state_component} ⊗ {b_state_component}"
    
    print("The initial entangled state is |\u03A8> = (1/\u221A5) \u03A3 |i>_A ⊗ |i>_B")
    print("After applying the gate Q to ququint A, the system is in a new superposition.")
    print("We measure ququint A. The probabilities of the outcomes |0> to |4> are:")
    # The normalization factor for the total state is 1/5.
    actual_probs = probs * (1/np.sum(probs))
    for i, p in enumerate(actual_probs):
        print(f"  P(A={i}) = {p:.2f}")
    
    print(f"\nThe most probable outcome is measuring ququint A in state |{most_probable_index}>.")
    print("The final state of the system after collapsing to this most probable outcome is:")
    print(final_equation)

solve_ququint_problem()