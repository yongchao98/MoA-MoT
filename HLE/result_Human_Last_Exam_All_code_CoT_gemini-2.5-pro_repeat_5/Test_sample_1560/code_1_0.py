import numpy as np

def solve_quantum_state():
    """
    Calculates and prints the possible final states and their probabilities
    after applying gate Q to the first of two entangled ququints and then
    measuring it.
    """
    # The gate Q is defined by its action on the basis states.
    # We can represent Q as a 5x5 matrix where the j-th column is Q|j>.
    # The problem defines:
    # Q|0> = (1/√2)(|1> + |2>) -> col 0: [0, 1, 1, 0, 0]
    # Q|1> = (1/√2)(|0> + |3>) -> col 1: [1, 0, 0, 1, 0]
    # Q|2> = (1/√2)(|1> + |4>) -> col 2: [0, 1, 0, 0, 1]
    # Q|3> = (1/√2)(|2> + |0>) -> col 3: [1, 0, 1, 0, 0]
    # Q|4> = (1/√2)(|3> + |2>) -> col 4: [0, 0, 1, 1, 0]
    # We define the integer part of the matrix M, and the 1/√2 factor separately.
    M = np.array([
        [0, 1, 0, 1, 0],
        [1, 0, 1, 0, 0],
        [1, 0, 0, 1, 1],
        [0, 1, 0, 0, 1],
        [0, 0, 1, 0, 0]
    ])
    
    # The matrix Q is M/√2
    Q = M / np.sqrt(2)

    print("After applying gate Q to the first ququint and measuring it, the possible outcomes are:")

    # The un-normalized state of B, correlated with A being in state |i>,
    # is derived from the i-th row of the Q matrix.
    # The probability of measuring |i>_A is (1/5) * ||row_i_of_Q||^2
    
    all_probs = []
    
    for i in range(5):
        # The i-th row of M gives the components of the un-normalized state of B.
        unnormalized_b_vector = M[i, :]
        
        # The squared norm is the sum of the squares of the components.
        norm_sq = np.sum(unnormalized_b_vector**2)
        
        # The probability of measuring |i>_A is (1/d) * ||row_i_of_Q||^2
        # which is (1/5) * (||row_i_of_M||^2 / 2) = norm_sq / 10
        probability = norm_sq / 10.0
        all_probs.append(probability)

        print(f"\n--- Outcome: ququint A is measured in state |{i}>_A ---")
        print(f"Probability of this outcome: {probability}")
        
        # Construct the string for the final state of B
        terms = []
        # Find the indices of non-zero elements
        ket_indices = np.where(unnormalized_b_vector != 0)[0]
        
        for k in ket_indices:
            coeff = int(unnormalized_b_vector[k])
            # Format the coefficient and ket state
            if norm_sq == 1:
                # Normalization factor is 1
                if coeff == 1:
                    terms.append(f"|{k}>")
                else:
                    terms.append(f"{coeff}|{k}>")
            else:
                # Normalization is 1/sqrt(norm_sq)
                if coeff == 1:
                    terms.append(f"(1/√{norm_sq})|{k}>")
                else: # Should not happen in this problem
                    terms.append(f"({coeff}/√{norm_sq})|{k}>")

        final_state_b_str = " + ".join(terms)
        print(f"Resulting state of ququint B: |ψ>_B = {final_state_b_str}")

    # Determine the most likely outcome for the final answer format
    most_likely_index = np.argmax(all_probs)
    # The final answer required is the single most likely measurement outcome.
    # print(f"\n<<<The most likely outcome is measuring ququint A in state |{most_likely_index}>, with probability {all_probs[most_likely_index]:.1f}.>>>")

solve_quantum_state()