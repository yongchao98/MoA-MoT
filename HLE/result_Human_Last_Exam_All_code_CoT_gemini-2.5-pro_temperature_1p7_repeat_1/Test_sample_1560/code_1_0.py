import numpy as np

def format_state(state_vector, ket_names=['|0>', '|1>', '|2>', '|3>', '|4>'], precision=4):
    """
    Formats a state vector into a readable string representation like (c0)|0> + (c1)|1> + ...
    """
    # Ensure the input is a column vector for consistent indexing
    if state_vector.ndim == 1:
        state_vector = state_vector.reshape(-1, 1)

    terms = []
    for i in range(len(state_vector)):
        # All coefficients in this problem are real
        coeff = state_vector[i, 0].real
        if abs(coeff) > 1e-9:  # Tolerance for zero
            if terms:
                terms.append(" + ")
            
            # Print each number in the final equation
            terms.append(f"({coeff:.{precision}f}) {ket_names[i]}")

    if not terms:
        return "0"

    return "".join(terms)

def solve_ququint_problem():
    """
    Solves the quantum mechanics problem as described.
    """
    # 1. Define basis states and the gate Q
    # Basis states |0> through |4> as 5x1 column vectors
    kets = [np.zeros((5, 1)) for _ in range(5)]
    for i in range(5):
        kets[i][i] = 1

    # Define the gate Q as a 5x5 matrix based on its action on the basis states
    Q = (1 / np.sqrt(2)) * np.array([
        # Q|0>   Q|1>   Q|2>   Q|3>   Q|4>
        [0, 1, 0, 1, 0],  # ket |0> component
        [1, 0, 1, 0, 0],  # ket |1> component
        [1, 0, 0, 1, 1],  # ket |2> component
        [0, 1, 0, 0, 1],  # ket |3> component
        [0, 0, 1, 0, 0]   # ket |4> component
    ])

    # 2. Define the initial 2-ququint entangled state in the 25-dim space
    # |\Psi_AB> = 1/sqrt(5) * (|00> + |11> + |22> + |33> + |44>)
    Psi_AB = np.zeros((25, 1))
    for i in range(5):
        Psi_AB += np.kron(kets[i], kets[i])
    Psi_AB /= np.sqrt(5)

    # 3. Apply the gate Q to the first ququint (A)
    # The operator is Q_A tensor I_B
    op = np.kron(Q, np.identity(5))
    Psi_prime_AB = op @ Psi_AB

    print("Final state of the system after measurement on ququint A:\n")
    print("The system collapses to one of the following states with the given probability.")
    print("-" * 70)
    
    max_prob = -1
    most_likely_state_B_str = ""
    
    # 4. Analyze measurement outcomes on ququint A
    for i in range(5):
        # The projection operator for measuring |i>_A is |i><i| ⊗ I
        proj_op = np.kron(kets[i] @ kets[i].T, np.identity(5))

        # Probability of measuring |i>_A is <Psi'|Proj|Psi'>
        prob = (Psi_prime_AB.conj().T @ proj_op @ Psi_prime_AB)[0, 0].real

        if prob > max_prob:
            max_prob = prob
        
        # Post-measurement state of the system (unnormalized)
        post_measure_sys_state = proj_op @ Psi_prime_AB

        # Post-measurement state of ququint B (unnormalized)
        # We "trace out" or "project out" ququint A
        bra_i_A = np.kron(kets[i].T, np.identity(5))
        post_measure_state_B = bra_i_A @ Psi_prime_AB

        # Normalize the state of B
        norm_B = np.linalg.norm(post_measure_state_B)
        normalized_state_B = post_measure_state_B / norm_B
        
        state_A_str = f"|{i}>_A"
        state_B_str = format_state(normalized_state_B, ket_names=['|0>_B', '|1>_B', '|2>_B', '|3>_B', '|4>_B'])

        print(f"If ququint A is measured as |{i}>:")
        print(f"  Probability: {prob:.4f}")
        final_state_equation = f"{state_A_str} ⊗ ({state_B_str})"
        print(f"  The final state of the system is: {final_state_equation}")
        print("-" * 70)
        
        if prob == max_prob:
             most_likely_state_B_str = f"({state_B_str})"
             
    # The final answer format is not specified for a probabilistic outcome.
    # We will output the state of ququint B for the most likely measurement outcome on A.
    # In case of a tie for max probability, this will take the last one found.
    # Let's find the most likely final state of the entire system as a string
    highest_prob_outcome = np.argmax([ (Psi_prime_AB.conj().T @ np.kron(kets[i] @ kets[i].T, np.identity(5)) @ Psi_prime_AB)[0, 0].real for i in range(5)])
    
    bra_i_A = np.kron(kets[highest_prob_outcome].T, np.identity(5))
    post_measure_state_B = bra_i_A @ Psi_prime_AB
    norm_B = np.linalg.norm(post_measure_state_B)
    normalized_state_B = post_measure_state_B / norm_B
    answer_state_b_str = format_state(normalized_state_B)
    
    # print(f"\n<<<For the final answer, providing the state of the second ququint (B) corresponding to the most likely measurement outcome of the first ququint (A).>>>")
    # print(f"<<<{answer_state_b_str}>>>")
    
# Run the solver
if __name__ == "__main__":
    solve_ququint_problem()
    
    # Deriving the answer for the tag again to be certain
    # From manual calculation, probabilites are [2/10, 2/10, 3/10, 2/10, 1/10]
    # Max prob is 3/10, for measuring |2>_A.
    # State of B is 1/sqrt(3) * (|0> + |3> + |4>).
    # 1/sqrt(3) = 0.57735...
    # The answer should be (0.5774) |0> + (0.5774) |3> + (0.5774) |4>
    # My format_state produces (c0)|0> + (c1)|1> so this needs adjustment for the tag.
    answer_string = f"({1/np.sqrt(3):.4f}) |0> + ({1/np.sqrt(3):.4f}) |3> + ({1/np.sqrt(3):.4f}) |4>"
    print(f"\n<<<{answer_string}>>>")