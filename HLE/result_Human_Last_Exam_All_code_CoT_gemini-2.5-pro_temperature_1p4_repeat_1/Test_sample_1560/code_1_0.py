import numpy as np

def final_equation_print(state_vector):
    """
    Formats a state vector into a symbolic string representation like
    (1/sqrt(N)) * (|i> + |j> + ...).
    This function also helps in fulfilling the requirement to output all numbers
    in the final equation by showing the normalization constant and the basis states involved.
    """
    # Find the indices of non-zero elements. Due to floating point arithmetic,
    # we check if a value is close to zero.
    non_zero_indices = np.where(np.isclose(state_vector, 0) == False)[0]
    
    num_terms = len(non_zero_indices)
    
    if num_terms == 0:
        return "0"
    
    # Create the sum of kets, e.g., "|1> + |3>"
    ket_terms = " + ".join([f"|{idx}>" for idx in non_zero_indices])

    if num_terms == 1:
        # The state is just a single basis state, e.g., "|2>"
        return ket_terms
    else:
        # The state is a superposition, e.g., "(1/sqrt(2)) * (|1> + |3>)"
        # The numbers in this equation are 1, 2, 1, and 3.
        return f"(1/sqrt({num_terms})) * ({ket_terms})"

def solve_ququint_problem():
    """
    Solves the main problem by applying gate Q to an entangled pair
    and calculating the post-measurement states.
    """
    # --- Step 1: Define the quantum system elements ---
    
    # Square root of 2 for normalization
    s_rt2 = np.sqrt(2)
    # Square root of 5 for normalization
    s_rt5 = np.sqrt(5)

    # Define the 5x5 matrix for the Q gate
    Q_gate = (1 / s_rt2) * np.array([
        [0, 1, 0, 1, 0],
        [1, 0, 1, 0, 0],
        [1, 0, 0, 1, 1],
        [0, 1, 0, 0, 1],
        [0, 0, 1, 0, 0]
    ], dtype=float)

    # Define the initial entangled state |Φ+> = 1/√5 * Σ |i>_A ⊗ |i>_B
    # This is a vector in a 25-dimensional Hilbert space (5x5).
    initial_state = np.zeros(25, dtype=float)
    for i in range(5):
        # The index in the 25D vector for state |i>⊗|i> is i*5 + i
        initial_state[i * 5 + i] = 1 / s_rt5

    # --- Step 2: Apply the gate to the first ququint ---
    
    # The operator for applying Q to the first ququint is Q ⊗ I
    identity_gate = np.identity(5)
    system_operator = np.kron(Q_gate, identity_gate)

    # Calculate the state of the system after the gate is applied
    state_after_gate = system_operator @ initial_state

    # --- Step 3: Calculate measurement outcomes ---

    # Reshape the 25x1 state vector into a 5x5 matrix of coefficients.
    # The element M[i, j] is the amplitude of the |i>_A ⊗ |j>_B state.
    coeff_matrix = state_after_gate.reshape((5, 5))

    print("After applying gate Q to the first ququint, the system is in a new state.")
    print("We now measure the first ququint. The possible outcomes and resulting system states are:\n")

    # Loop over the 5 possible measurement outcomes for ququint A (|0> through |4>)
    for i in range(5):
        # The i-th row of the matrix contains the amplitudes for ququint B
        # given that ququint A was measured in state |i>.
        state_b_unnormalized = coeff_matrix[i, :]
        
        # The probability of measuring |i>_A is the sum of the squared magnitudes
        # of the amplitudes in this row (the squared norm of the vector).
        prob_i = np.sum(np.abs(state_b_unnormalized)**2)
        
        # We only describe outcomes with a non-zero probability of occurring.
        if not np.isclose(prob_i, 0):
            # To find the state of ququint B, we normalize its state vector.
            norm_factor = np.sqrt(prob_i)
            state_b_normalized = state_b_unnormalized / norm_factor
            
            # --- Step 4: Print the results ---
            print(f"--- Outcome: Ququint A is measured as |{i}> ---")
            # Probabilities are fractions with denominator 10.
            prob_numerator = int(round(prob_i * 10))
            print(f"  Probability: {prob_i:.2f} ({prob_numerator}/10)")
            
            # The final state of the entire system is |i>_A ⊗ |ψ>_B
            print(f"  The final system state is |{i}>_A ⊗ |ψ>_B, where:")
            
            # Print the equation for the state of ququint B
            final_equation = final_equation_print(state_b_normalized)
            print(f"  |ψ>_B = {final_equation}\n")

# Run the simulation and print the results
solve_ququint_problem()