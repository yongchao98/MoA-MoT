import numpy as np

def format_ket(vector, system_label, precision=1e-9):
    """Formats a state vector into a readable ket notation string."""
    terms = []
    
    # Find unique non-zero coefficients to group terms
    coeffs = np.unique(vector[np.abs(vector) > precision])
    
    if len(coeffs) == 0:
        return "0"

    # Handle the case where the state is a single basis vector
    if len(coeffs) == 1 and np.isclose(coeffs[0], 1.0):
        idx = np.where(np.isclose(vector, 1.0))[0][0]
        return f"|{idx}>_{system_label}"

    for c in coeffs:
        indices = np.where(np.isclose(vector, c))[0]
        
        # Format the coefficient
        coeff_str = ""
        # Check for common fractions like 1/sqrt(N)
        if np.isclose(c**2, 1/2):
            coeff_str = f"(1/sqrt(2))"
        elif np.isclose(c**2, 1/3):
            coeff_str = f"(1/sqrt(3))"
        elif np.isclose(c, 1.0):
            coeff_str = ""
        elif np.isclose(c, -1.0):
            coeff_str = "-"
        else:
            coeff_str = f"{c:.3f}"

        # Format the ket part
        if len(indices) == 1:
            ket_str = f"|{indices[0]}>_{system_label}"
        else:
            ket_sum = " + ".join([f"|{i}>_{system_label}" for i in indices])
            ket_str = f"({ket_sum})"
        
        if coeff_str:
            terms.append(f"{coeff_str}{ket_str}")
        else:
            terms.append(ket_str)
            
    return " + ".join(terms).replace("+ -", "- ")

def solve_ququint_measurement():
    """
    Calculates and prints the final states of the two-ququint system after measurement.
    """
    d = 5  # Dimension of the ququint system

    # Define the gate Q as a matrix
    # Q|0> = 1/sqrt(2) * (|1> + |2>) -> col 0
    # Q|1> = 1/sqrt(2) * (|0> + |3>) -> col 1
    # etc.
    Q = (1 / np.sqrt(2)) * np.array([
        [0, 1, 0, 1, 0],
        [1, 0, 1, 0, 0],
        [1, 0, 0, 1, 1],
        [0, 1, 0, 0, 1],
        [0, 0, 1, 0, 0]
    ], dtype=float)

    # The state after applying Q to ququint A is
    # |Psi'> = (1/sqrt(5)) * sum_{i=0 to 4} (Q|i>_A) x |i>_B
    # We can find the unnormalized state of B for each state of A
    # by collecting terms.
    
    # phi_b_unnormalized[j] will store the unnormalized state of B
    # when A is in state |j>
    phi_b_components = [np.zeros(d) for _ in range(d)]

    for i in range(d):
        # This is Q|i>
        q_ket_i = Q[:, i]
        for j in range(d):
            # The coefficient of |j>_A in Q|i>_A is q_ket_i[j]
            # This gets paired with |i>_B
            phi_b_components[j] += q_ket_i[j] * np.eye(d)[i]

    print("The system starts in the entangled state:")
    print("|Psi_entangled> = (1/sqrt(5)) * (|0>_A|0>_B + |1>_A|1>_B + |2>_A|2>_B + |3>_A|3>_B + |4>_A|4>_B)\n")
    
    print("After applying gate Q to ququint A, the state is rearranged for measurement:")
    # The overall normalization factor is (1/sqrt(5))
    # So we calculate probabilities and states from here.
    
    print("\nUpon measuring ququint A, the system collapses to one of the following outcomes:\n")

    total_prob = 0
    for j in range(d):
        # The full unnormalized state of B is (1/sqrt(5)) * phi_b_components[j]
        # The probability is the norm squared of this vector.
        prob = (1/5) * np.linalg.norm(phi_b_components[j])**2
        total_prob += prob
        
        # The final state of B is the normalized version of phi_b_components[j]
        norm_phi_b = np.linalg.norm(phi_b_components[j])
        if norm_phi_b > 1e-9:
            final_state_b = phi_b_components[j] / norm_phi_b
        else:
            final_state_b = np.zeros(d)

        print(f"--- Outcome |{j}>_A ---")
        print(f"Probability = {prob:.3f} (or {int(prob*10)}/10)")
        
        final_state_b_str = format_ket(final_state_b, "B")
        print(f"Final State = |{j}>_A ⊗ {final_state_b_str}\n")

    print(f"Total probability check: {total_prob:.3f}")

solve_ququint_measurement()