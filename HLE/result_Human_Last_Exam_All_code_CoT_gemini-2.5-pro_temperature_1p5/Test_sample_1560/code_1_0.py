import numpy as np
import math

def solve_ququint_state():
    """
    Calculates the final state of a ququint system after applying a given quantum gate.

    The steps are:
    1. Define the initial state vector and the gate matrix from the problem description.
    2. Apply the gate to the state via matrix-vector multiplication.
    3. Calculate the norm of the resulting unnormalized state.
    4. Construct and print the string representation of the final, normalized state.
    """
    # We define the vector of coefficients for the initial state |ψ⟩.
    # |ψ⟩ = 1/sqrt(5) * [1, 1, 1, 1, 1]^T. We can work with the integer part for now.
    psi_coeffs = np.array([1, 1, 1, 1, 1])

    # We define the matrix M based on the definition of the gate Q,
    # where Q = (1/sqrt(2)) * M.
    M = np.array([
        [0, 1, 0, 1, 0],
        [1, 0, 1, 0, 0],
        [1, 0, 0, 1, 1],
        [0, 1, 0, 0, 1],
        [0, 0, 1, 0, 0]
    ])

    # Apply the integer part of the gate to the integer part of the state vector.
    # The overall normalization constant will be handled at the end.
    final_coeffs_unnormalized = M @ psi_coeffs

    # Calculate the squared norm of the resulting coefficient vector. This will
    # be used in the denominator of the final coefficients.
    norm_squared = np.sum(final_coeffs_unnormalized**2)
    norm_val = int(norm_squared)

    # Build the string representation of the final state equation.
    # The final state is Σ (c_i / sqrt(norm_val)) |i⟩
    # where c_i are the components of final_coeffs_unnormalized.
    output_string = "The final state is: "
    is_first_term = True
    for i, coeff in enumerate(final_coeffs_unnormalized):
        # We can skip terms with a zero coefficient
        if coeff == 0:
            continue
        
        if not is_first_term:
            output_string += " + "
        
        # Add the term for the current basis state to the equation string
        output_string += f"({int(coeff)}/sqrt({norm_val}))|{i}⟩"
        is_first_term = False

    print(output_string)

solve_ququint_state()