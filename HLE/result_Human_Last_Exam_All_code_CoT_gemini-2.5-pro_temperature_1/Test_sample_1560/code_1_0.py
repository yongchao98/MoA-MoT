import numpy as np

def solve_ququint_problem():
    """
    Calculates the final state of a ququint system after applying a custom gate Q.
    """
    # Define the square root of 2 for convenience
    sqrt2 = np.sqrt(2)

    # The gate Q is defined by its action on the basis states.
    # We can represent Q as a matrix where the columns are the results of Q|i>.
    # Q|0> = 1/sqrt(2) * (|1> + |2>) -> [0, 1, 1, 0, 0] / sqrt(2)
    # Q|1> = 1/sqrt(2) * (|0> + |3>) -> [1, 0, 0, 1, 0] / sqrt(2)
    # Q|2> = 1/sqrt(2) * (|1> + |4>) -> [0, 1, 0, 0, 1] / sqrt(2)
    # Q|3> = 1/sqrt(2) * (|2> + |0>) -> [1, 0, 1, 0, 0] / sqrt(2)
    # Q|4> = 1/sqrt(2) * (|3> + |2>) -> [0, 0, 1, 1, 0] / sqrt(2)

    Q_matrix = (1 / sqrt2) * np.array([
        [0, 1, 0, 1, 0],
        [1, 0, 1, 0, 0],
        [1, 0, 0, 1, 1],
        [0, 1, 0, 0, 1],
        [0, 0, 1, 0, 0]
    ])

    # The initial state |ψ> = 1/sqrt(5) * (|0> + |1> + |2> + |3> + |4>)
    psi_initial_vector = (1 / np.sqrt(5)) * np.ones(5)

    # Calculate the final state |ψ'> by applying the gate Q to |ψ>
    # |ψ'> = Q |ψ>
    psi_final_vector = Q_matrix @ psi_initial_vector

    # The resulting vector has a common factor of 1/sqrt(10).
    # Let's find the integer numerators for the coefficients.
    # psi_final_vector = (1/sqrt(10)) * [c0, c1, c2, c3, c4]
    # So, [c0, c1, c2, c3, c4] = psi_final_vector * sqrt(10)
    numerators = np.round(psi_final_vector * np.sqrt(10)).astype(int)
    denominator_str = "sqrt(10)"

    # Construct the output string for the final state equation
    final_state_parts = []
    for i, num in enumerate(numerators):
        if num != 0:
            # Handle the case where the numerator is 1
            if num == 1:
                coeff_str = f"1/{denominator_str}"
            else:
                coeff_str = f"{num}/{denominator_str}"
            final_state_parts.append(f"({coeff_str})|{i}>")

    final_state_equation = " + ".join(final_state_parts)

    print("The final state |ψ'> before measurement is:")
    print(f"|ψ'> = {final_state_equation}")

    # According to the problem's measurement rule, the probability of collapsing
    # to state |i> is the square of the amplitude |a_i|^2.
    # Let's calculate the sum of these "probabilities".
    # A standard unitary gate would preserve the norm, and the sum would be 1.
    sum_of_probabilities = np.sum(np.abs(psi_final_vector)**2)
    
    print("\nFor context, the sum of the squared amplitudes (the 'total probability' under the given rules) is:")
    print(sum_of_probabilities)


solve_ququint_problem()
<<<2.2>>>