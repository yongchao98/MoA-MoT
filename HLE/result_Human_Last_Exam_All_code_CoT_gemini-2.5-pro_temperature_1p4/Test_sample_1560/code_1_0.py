import numpy as np
import math

def solve_ququint_state():
    """
    Calculates the final state of the ququint system and prints the result.
    """
    # 1. Define the gate Q as a matrix.
    # The action of Q on the basis vectors gives the columns of the matrix.
    # For example, Q|0> = (1/sqrt(2))(|1> + |2>) corresponds to the column [0, 1, 1, 0, 0]^T.
    # The factor 1/sqrt(2) is common to all definitions.
    q_matrix = (1/math.sqrt(2)) * np.array([
        [0, 1, 0, 1, 0],
        [1, 0, 1, 0, 0],
        [1, 0, 0, 1, 1],
        [0, 1, 0, 0, 1],
        [0, 0, 1, 0, 0]
    ])

    # 2. Define the initial state |psi> as a vector.
    # |psi> = 1/sqrt(5) * (|0> + |1> + |2> + |3> + |4>)
    psi_vector = (1/math.sqrt(5)) * np.array([1, 1, 1, 1, 1])

    # 3. Calculate the final state |psi'> by applying Q to |psi>.
    # |psi'> = Q |psi>
    final_state_vector = q_matrix @ psi_vector

    # The resulting vector's coefficients are with respect to the basis |0>, |1>, ...
    # The overall normalization factor is (1/sqrt(2)) * (1/sqrt(5)) = 1/sqrt(10).
    # We can get the integer coefficients by multiplying by sqrt(10).
    unnormalized_coeffs = final_state_vector * math.sqrt(10)
    
    # We round to handle potential floating point inaccuracies and convert to int
    coefficients = [int(round(c)) for c in unnormalized_coeffs]

    # 4. Format and print the final state equation.
    # The prompt requires outputting each number in the final equation.
    print("The final state of the system before measurement is described by the equation:")
    
    equation_parts = []
    for i, coeff in enumerate(coefficients):
        if coeff != 0:
            equation_parts.append(f"{coeff}/sqrt(10) |{i}>")

    final_equation = "Final State = " + " + ".join(equation_parts)
    print(final_equation)
    
    # As per the plan, we also calculate the expectation value for the final answer.
    # Unnormalized probabilities are the square of the coefficients.
    unnormalized_probs = np.array(coefficients)**2
    # The normalization constant is the sum of these squared coefficients.
    norm_sum = np.sum(unnormalized_probs)
    # Normalized probabilities
    probabilities = unnormalized_probs / norm_sum
    # Measurement outcomes
    outcomes = np.arange(len(coefficients))
    # Expectation value
    expectation_value = np.sum(outcomes * probabilities)
    
    # This value will be provided in the final answer block.
    # print(f"\nThe expectation value of a measurement is: {expectation_value:.3f}")


solve_ququint_state()