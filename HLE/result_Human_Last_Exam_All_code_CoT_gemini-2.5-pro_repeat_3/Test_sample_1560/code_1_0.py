import numpy as np
import math

def solve_ququint_state():
    """
    Calculates the final state of a ququint system after a gate operation.
    """
    # 1. Define the matrix for the gate Q (without the 1/sqrt(2) factor).
    # The columns of the matrix are the results of Q acting on the basis states |0> through |4>.
    Q_matrix = np.array([
        [0, 1, 0, 1, 0],  # Q|0> is col 0, Q|1> is col 1, etc.
        [1, 0, 1, 0, 0],
        [1, 0, 0, 1, 1],
        [0, 1, 0, 0, 1],
        [0, 0, 1, 0, 0]
    ])

    # 2. Define the initial state |psi> (without the 1/sqrt(5) factor).
    psi_vector = np.array([1, 1, 1, 1, 1])

    # 3. Apply the gate's matrix to the state's vector.
    # Note: The initial normalization factors (1/sqrt(2) and 1/sqrt(5)) are carried through,
    # but it's simpler to apply the matrices without them and normalize at the end.
    final_state_unnormalized = Q_matrix @ psi_vector

    # 4. Calculate the squared norm of the resulting unnormalized vector.
    # The norm is the square root of the sum of the squares of the coefficients.
    norm_squared = np.sum(final_state_unnormalized**2)
    
    # 5. Construct the string for the final equation.
    # The final state is |psi'> = (1/sqrt(norm_squared)) * (c0|0> + c1|1> + ...)
    # where c_i are the components of the unnormalized final state vector.
    
    coeffs = final_state_unnormalized
    
    # Build the part of the equation inside the parentheses
    terms = []
    for i in range(len(coeffs)):
        # Add a plus sign for non-negative coefficients after the first term
        sign = " + " if coeffs[i] >= 0 and i > 0 else " "
        term_str = f"{sign}{coeffs[i]}|{i}>"
        terms.append(term_str)
    
    # Assemble the final equation string
    equation = f"|psi'> = (1/sqrt({norm_squared})) * ({''.join(terms).lstrip(' + ')})"

    print("The final state of the system before measurement is described by the equation:")
    print(equation)

# Execute the function
solve_ququint_state()