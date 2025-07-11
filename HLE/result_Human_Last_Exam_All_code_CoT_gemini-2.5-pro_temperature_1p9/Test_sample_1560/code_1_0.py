import numpy as np

def solve_quantum_state():
    """
    Calculates the final state of a ququint system after a gate operation.
    """
    # 1. Define the gate Q as a 5x5 matrix.
    # The definitions of Q|i> give the columns of the matrix.
    # For example, Q|0> = (1/sqrt(2))(|1> + |2>) corresponds to the column [0, 1, 1, 0, 0]^T / sqrt(2).
    q_matrix = (1 / np.sqrt(2)) * np.array([
        [0, 1, 0, 1, 0],
        [1, 0, 1, 0, 0],
        [1, 0, 0, 1, 1],
        [0, 1, 0, 0, 1],
        [0, 0, 1, 0, 0]
    ])

    # 2. Define the initial state |psi> as a 5-element vector.
    # |psi> = 1/sqrt(5) * (|0> + |1> + |2> + |3> + |4>)
    psi_vector = (1 / np.sqrt(5)) * np.ones(5)

    # 3. Apply the gate Q to the state |psi> by matrix-vector multiplication.
    final_state_vector = np.dot(q_matrix, psi_vector)

    # 4. Print the final state in the standard ket notation.
    # The problem asks to output each number in the final equation.
    print("The final state |psi'> = Q|psi> is:")
    
    equation_parts = []
    for i in range(len(final_state_vector)):
        # Format the coefficient to four decimal places for clarity.
        coeff = final_state_vector[i]
        equation_parts.append(f"({coeff:.4f})|{i}>")
    
    final_equation = " + ".join(equation_parts)
    print(final_equation)

    # Return the string for the final answer block
    return final_equation

# Execute the function
final_state_equation = solve_quantum_state()

# The final answer in the requested format
# print(f"<<<{final_state_equation}>>>")