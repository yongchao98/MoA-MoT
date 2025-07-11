import numpy as np

def solve_switch_problem():
    """
    Calculates the expected number of rounds for the switch system to return to its initial state.
    """
    # Define influence sets based on the problem (1-based indexing for people)
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

    # Create the influence matrix A (using 0-based indexing for matrix)
    # A[j][i] = 1 if person i+1 influences person j+1
    A = np.zeros((8, 8), dtype=int)
    for i in range(1, 9):
        for j in influence_sets[i]:
            A[j - 1, i - 1] = 1

    # The transition matrix T = I + A (all operations are modulo 2)
    T = (np.identity(8, dtype=int) + A) % 2

    num_states = 2**8
    # Dictionary to store the cycle length for each state, to avoid re-computation.
    cycle_lengths = {}

    for i in range(num_states):
        # Create the initial state vector s0 as a tuple from the integer i.
        # A tuple is hashable and can be used as a dictionary key.
        s0 = tuple(int(x) for x in np.binary_repr(i, width=8))

        if s0 in cycle_lengths:
            continue

        # Simulate the rounds to find the cycle length
        path = [s0]
        s_current = tuple(T.dot(np.array(s0)) % 2)
        while s_current != s0:
            path.append(s_current)
            s_current = tuple(T.dot(np.array(s_current)) % 2)
        
        # The length of the path is the cycle length.
        length = len(path)
        
        # All states in the discovered cycle have the same return time.
        for state_in_cycle in path:
            cycle_lengths[state_in_cycle] = length

    # Count the number of states for each cycle length
    counts = {}
    for length in cycle_lengths.values():
        counts[length] = counts.get(length, 0) + 1
    
    # Sort counts by key for ordered printing
    sorted_counts = sorted(counts.items())

    # Calculate the sum for the expectation formula
    total_sum_of_rounds = sum(cycle_lengths.values())
    expected_value = total_sum_of_rounds / num_states
    
    print("The number of states for each cycle length:")
    term_strings = []
    value_strings = []
    for length, count in sorted_counts:
        print(f"Length {length}: {count}")
        term_strings.append(f"{length} * {count}")
        value_strings.append(str(length * count))

    print("\nCalculating the expected value E[R]:")
    
    # Printing the detailed equation as requested
    equation_part1 = f"({ ' + '.join(term_strings) }) / {num_states}"
    equation_part2 = f"({ ' + '.join(value_strings) }) / {num_states}"
    equation_part3 = f"{total_sum_of_rounds} / {num_states}"

    print(f"E[R] = {equation_part1}")
    print(f"     = {equation_part2}")
    print(f"     = {equation_part3}")
    print(f"     = {expected_value:.2f}")


solve_switch_problem()
<<<7.46>>>