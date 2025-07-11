import numpy as np

def state_to_string(state_vector, system_label=""):
    """
    Converts a state vector into a human-readable string format,
    like "0.707|1> + 0.707|3>".
    """
    # A small number to check if a coefficient is effectively zero
    epsilon = 1e-9
    components = []
    for i, amp in enumerate(state_vector.flatten()):
        # Only include terms with non-zero amplitudes
        if abs(amp.real) > epsilon or abs(amp.imag) > epsilon:
            # Format the coefficient
            if abs(amp.imag) < epsilon: # Real number
                coeff = f"{amp.real:.3f}"
            else: # Complex number
                coeff = f"({amp.real:.3f}{amp.imag:+.3f}j)"
            
            # Simplify for coefficients that are 1.0
            if coeff == "1.000":
                components.append(f"|{i}>{system_label}")
            else:
                components.append(f"{coeff}|{i}>{system_label}")
    
    if not components:
        return "0"
        
    # Join the components with ' + ' and clean up formatting
    return " + ".join(components).replace("+ -", "- ")

def solve_ququint_measurement():
    """
    Solves the ququint measurement problem as described.
    """
    # 1. Define the basis states and the dimension of the system
    d = 5
    basis_states = [np.zeros((d, 1)) for _ in range(d)]
    for i in range(d):
        basis_states[i][i] = 1

    # 2. Define the gate Q as a matrix
    # The columns of the matrix are Q|0>, Q|1>, ...
    sqrt2_inv = 1 / np.sqrt(2)
    Q0 = sqrt2_inv * (basis_states[1] + basis_states[2])
    Q1 = sqrt2_inv * (basis_states[0] + basis_states[3])
    Q2 = sqrt2_inv * (basis_states[1] + basis_states[4])
    Q3 = sqrt2_inv * (basis_states[2] + basis_states[0])
    Q4 = sqrt2_inv * (basis_states[3] + basis_states[2])
    
    Q = np.hstack([Q0, Q1, Q2, Q3, Q4])

    # 3. Define the initial entangled state
    # |Psi_entangled> = 1/sqrt(5) * sum(|i>_A |i>_B)
    sqrt5_inv = 1 / np.sqrt(5)
    psi_entangled = np.zeros((d*d, 1))
    for i in range(d):
        # The state |i>_A |i>_B is the Kronecker product of |i> and |i>
        ket_i_A = basis_states[i]
        ket_i_B = basis_states[i]
        psi_entangled += np.kron(ket_i_A, ket_i_B)
    psi_entangled *= sqrt5_inv

    # 4. Apply the gate Q to the first ququint (A)
    # The operator is Q_A tensor I_B
    identity_B = np.identity(d)
    op = np.kron(Q, identity_B)
    psi_new = op @ psi_entangled

    print("After applying gate Q to the first ququint, we analyze the measurement outcomes for ququint A.\n")

    # 5. Analyze the measurement outcomes for ququint A
    total_prob = 0
    for i in range(d):
        # Extract the part of the state corresponding to measuring |i>_A
        # This is a vector of coefficients for ququint B's state.
        # The full state vector psi_new is ordered as c_00, c_01.. c_04, c_10, c_11..
        # so for |i>_A, B's state is in the slice [i*d : (i+1)*d]
        state_B_unnormalized = psi_new[i*d : (i+1)*d]
        
        # The probability of this outcome is the squared norm of this sub-vector
        prob = np.linalg.norm(state_B_unnormalized)**2
        total_prob += prob

        # The state of ququint B collapses to the normalized version of this vector
        if prob > 1e-9:
            state_B_normalized = state_B_unnormalized / np.sqrt(prob)
        else:
            state_B_normalized = np.zeros((d, 1))

        # 6. Print the results for this outcome
        # Using fractions for probabilities where they are neat (multiples of 1/10)
        prob_frac_numerator = int(round(prob * 10))
        prob_frac = f"{prob_frac_numerator}/10"

        print(f"Outcome {i+1}: Measure ququint A in state |{i}>")
        print(f"Probability: {prob:.2f} (or {prob_frac})")
        
        # Format the final state equation
        final_state_B_str = state_to_string(state_B_normalized, "_B")
        print(f"The final system state is: |{i}>_A ⊗ ( {final_state_B_str} )")
        print("-" * 50)
    
    print(f"Total probability check: {total_prob:.2f}")


# Run the simulation and print the final answer
solve_ququint_measurement()