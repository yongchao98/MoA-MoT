import numpy as np

def solve_ququint_state():
    """
    Calculates the final state of a ququint system after applying gate Q.
    """
    # Define the gate Q as a 5x5 matrix.
    # The columns of the matrix are the results of Q acting on the basis states |0> through |4>.
    # Q|0> = 1/sqrt(2) * (|1> + |2>) -> [0, 1, 1, 0, 0]
    # Q|1> = 1/sqrt(2) * (|0> + |3>) -> [1, 0, 0, 1, 0]
    # Q|2> = 1/sqrt(2) * (|1> + |4>) -> [0, 1, 0, 0, 1]
    # Q|3> = 1/sqrt(2) * (|2> + |0>) -> [1, 0, 1, 0, 0]
    # Q|4> = 1/sqrt(2) * (|3> + |2>) -> [0, 0, 1, 1, 0]
    
    Q_matrix = (1 / np.sqrt(2)) * np.array([
        [0, 1, 0, 1, 0],
        [1, 0, 1, 0, 0],
        [1, 0, 0, 1, 1],
        [0, 1, 0, 0, 1],
        [0, 0, 1, 0, 0]
    ])

    # Define the initial state psi = 1/sqrt(5) * (|0> + |1> + |2> + |3> + |4>)
    psi_vector = (1 / np.sqrt(5)) * np.ones(5)

    # Calculate the final state by multiplying the Q matrix and the psi vector.
    psi_final_vector = Q_matrix @ psi_vector
    
    # Construct the final equation string.
    # The format will be: c_0|0> + c_1|1> + ...
    # where c_i are the coefficients from the final vector.
    
    equation_parts = []
    for i, coeff in enumerate(psi_final_vector):
        # We present the coefficients in their symbolic form, derived from the calculation.
        # The vector before normalization is [2, 2, 3, 2, 1] / sqrt(10).
        numerator = int(round(coeff * np.sqrt(10)))
        if numerator != 0:
            term = f"({numerator}/sqrt(10))|{i}>"
            equation_parts.append(term)
    
    final_equation = " + ".join(equation_parts)

    print("The final state of the system before measurement is:")
    print(final_equation)

    # To be more explicit with the numbers as requested.
    print("\nWhich numerically is:")
    numeric_equation_parts = []
    for i, coeff in enumerate(psi_final_vector):
         if abs(coeff) > 1e-9: # a small tolerance for zero
            numeric_equation_parts.append(f"({coeff:.4f})|{i}>")
    numeric_equation = " + ".join(numeric_equation_parts)
    print(numeric_equation)

# Run the solver
solve_ququint_state()