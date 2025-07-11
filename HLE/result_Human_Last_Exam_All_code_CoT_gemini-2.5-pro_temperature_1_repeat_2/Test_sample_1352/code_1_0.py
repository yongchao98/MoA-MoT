def solve_switch_problem():
    """
    Calculates the expected number of rounds for the switch system to return to its initial state.
    """
    # Define the influence sets as given in the problem description.
    # Person numbers are 1-based.
    influence_sets = {
        1: {2, 4, 6, 7},
        2: {3, 5, 6, 8},
        3: {4, 6},
        4: {5},
        5: {6, 8},
        6: {7},
        7: {8},
        8: {}
    }

    # Create the 8x8 transformation matrix M for operations over F_2.
    # Matrix indices are 0-based (0 to 7).
    M = [[0] * 8 for _ in range(8)]

    # The new state s_out is M * s_in.
    # s_out[j] = s_in[j] + sum over k where j is in I_k (s_in[k]).
    # This means M[j][j] = 1, and M[j][k] = 1 if j is influenced by k.
    for i in range(8):
        M[i][i] = 1  # Diagonal elements

    for k_person, influenced_people in influence_sets.items():
        for j_person in influenced_people:
            # Row index j corresponds to person j+1
            # Column index k corresponds to person k+1
            M[j_person - 1][k_person - 1] = 1

    def mat_vec_mul_mod2(matrix, vector):
        """Helper function for matrix-vector multiplication modulo 2."""
        n = len(vector)
        result = [0] * n
        for i in range(n):
            sum_val = 0
            for j in range(n):
                sum_val += matrix[i][j] * vector[j]
            result[i] = sum_val % 2
        return tuple(result)

    # Iterate through all 2^8 = 256 possible initial states.
    num_states = 256
    total_rounds_sum = 0

    for i in range(num_states):
        # Generate the initial state S_0 as an 8-bit binary vector.
        initial_state = tuple(int(bit) for bit in format(i, '08b'))

        # Find the cycle length for this state.
        current_state = initial_state
        rounds = 0
        while True:
            # Apply one round of transformation.
            current_state = mat_vec_mul_mod2(M, list(current_state))
            rounds += 1
            if current_state == initial_state:
                total_rounds_sum += rounds
                break
    
    # Calculate the expected value E[R].
    expected_R = total_rounds_sum / num_states

    # Print the numbers used in the final equation and the result.
    print(f"The sum of cycle lengths for all 256 states is: {total_rounds_sum}")
    print(f"The total number of states is: {num_states}")
    print("The final equation for the expected value is:")
    print(f"{total_rounds_sum} / {num_states} = {expected_R}")
    print(f"The expected value E[R] rounded to 2 decimal places is: {expected_R:.2f}")

solve_switch_problem()
<<<4.52>>>