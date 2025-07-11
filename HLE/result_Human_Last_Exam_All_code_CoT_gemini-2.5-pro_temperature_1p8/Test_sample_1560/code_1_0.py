import numpy as np

def solve_quantum_system():
    """
    Solves the ququint quantum system problem.
    1. Defines the gate Q and the initial entangled state.
    2. Applies the gate Q to the first ququint.
    3. Analyzes the resulting state to determine the outcomes of a measurement
       on the first ququint.
    4. Prints the final results as readable equations.
    """
    # Dimension of the ququint system
    d = 5

    # Define basis states |0> to |4> as column vectors
    basis = [np.zeros((d, 1)) for _ in range(d)]
    for i in range(d):
        basis[i][i] = 1

    # 1. Construct the matrix for gate Q based on its definition
    # Q|0⟩ = 1/√2 * (|1⟩ + |2⟩)
    # Q|1⟩ = 1/√2 * (|0⟩ + |3⟩)
    # Q|2⟩ = 1/√2 * (|1⟩ + |4⟩)
    # Q|3⟩ = 1/√2 * (|2⟩ + |0⟩)
    # Q|4⟩ = 1/√2 * (|3⟩ + |2⟩)
    Q_matrix = (1 / np.sqrt(2)) * np.array([
        [0, 1, 0, 1, 0],
        [1, 0, 1, 0, 0],
        [1, 0, 0, 1, 1],
        [0, 1, 0, 0, 1],
        [0, 0, 1, 0, 0]
    ])

    # 2. Define the initial entangled state
    # |Ψ⟩ = 1/√5 * Σ (|i⟩_A ⊗ |i⟩_B)
    psi_entangled = np.zeros((d * d, 1))
    for i in range(d):
        psi_entangled += np.kron(basis[i], basis[i])
    psi_entangled *= (1 / np.sqrt(d))

    # 3. Apply gate Q to the first ququint (A).
    # The operator is (Q ⊗ I), where I is the identity matrix.
    op = np.kron(Q_matrix, np.identity(d))
    final_state_vec = op @ psi_entangled

    # 4. Analyze and print the state BEFORE measurement
    # This state determines all possible measurement outcomes.
    # The final state vector has an overall normalization factor of 1/√(2*5) = 1/√10
    print("The state of the system after applying gate Q to ququint A (before measurement) is:")
    print("|Ψ'⟩ = (1 / sqrt(10)) * [")

    # Group terms by ququint A's basis states to make it readable
    final_state_terms = {}
    for i in range(d):
        # For each |i>_A, find the corresponding state of B
        state_B_parts = []
        for j in range(d):
            # Coefficient for the |i⟩|j⟩ term
            coeff = final_state_vec[i * d + j][0]
            if not np.isclose(coeff, 0):
                state_B_parts.append(f"|{j}⟩_B")
        if state_B_parts:
            final_state_terms[i] = state_B_parts

    output_lines = []
    for i in sorted(final_state_terms.keys()):
        b_states_str = " + ".join(final_state_terms[i])
        if len(final_state_terms[i]) > 1:
            b_states_str = f"({b_states_str})"
        output_lines.append(f"    |{i}⟩_A ⊗ {b_states_str}")

    print("   + \n".join(output_lines))
    print("]\n")

    # 5. Determine the final state AFTER measurement by analyzing outcomes
    print("A measurement on the first ququint (A) will collapse the system into one of the following states:")
    
    total_prob = 0
    for i in sorted(final_state_terms.keys()):
        num_terms = len(final_state_terms[i])
        prob = num_terms / 10.0
        total_prob += prob

        b_states_str = " + ".join(final_state_terms[i])

        # Create the normalized state for ququint B
        if num_terms > 1:
            normalized_b_state = f"(1/sqrt({num_terms})) * ({b_states_str})"
        else:
            normalized_b_state = b_states_str # Already normalized

        # Print the outcome, probability, and resulting final state
        print("---")
        print(f"Outcome: Ququint A is measured in state |{i}⟩_A")
        print(f"Probability: {num_terms}/10 = {prob:.1f}")
        print(f"Resulting Final State: |{i}⟩_A ⊗ {normalized_b_state}")

    # Final check that probabilities sum to 1
    print("---")
    print(f"Total Probability Check: {total_prob:.1f}")

if __name__ == '__main__':
    solve_quantum_system()
    # The task implies that the final state is the result of a measurement.
    # Because measurement is probabilistic, the complete description involves
    # all possible collapsed states and their probabilities.
    # The most comprehensive answer is the detailed breakdown printed above.
    # As there is no single deterministic "final state", a specific value like a number or a letter
    # cannot represent the solution. The "answer" is the full description of outcomes.
    # We choose the probability of the most likely outcome as the representative numerical answer.
    # P(A=2) = 3/10 = 0.3
    print("\n<<<0.3>>>")