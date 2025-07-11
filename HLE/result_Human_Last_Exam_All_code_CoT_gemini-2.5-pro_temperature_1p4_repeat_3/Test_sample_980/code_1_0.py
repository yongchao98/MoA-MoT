import itertools

def find_next_number(seq):
    """
    Analyzes the sequence by grouping identical consecutive numbers
    and extrapolating the pattern to find the next number.
    """
    # Step 1 & 2: Group the sequence and create a list of (value, count) pairs.
    groups = [(key, len(list(group))) for key, group in itertools.groupby(seq)]

    print("Original sequence:", seq)
    print("Groups of (value, count):", groups)
    print("-" * 20)
    
    # Step 3 & 4: Analyze the pattern at the end of the sequence.
    # The pattern of the last four values is (36, 124, 36, 124).
    last_four_values = [g[0] for g in groups[-4:]]
    
    # The logic is based on pattern continuation.
    # Pattern: ..., value_A, value_B, value_A, value_B, ...
    # The next value is predicted to be value_A.
    predicted_next_value = last_four_values[-2]

    # The final step is to formulate the "equation" as requested.
    # The equation is a statement of the pattern logic.
    val1 = last_four_values[0]
    val2 = last_four_values[1]
    val3 = last_four_values[2]
    val4 = last_four_values[3]
    
    # Step 5: Print the reasoning and the final answer.
    print("The last four unique values in the sequence are {}, {}, {}, {}.".format(val1, val2, val3, val4))
    print("This reveals an alternating pattern.")
    print("The equation to find the next number is based on continuing this pattern:")
    print("{} (value from 2 steps back) -> {} (predicted next value)".format(val3, predicted_next_value))
    print("\nTherefore, the next number in the sequence is:")
    print(predicted_next_value)


# The given sequence
sequence = [
    111, 142, 111, 41, 67, 67, 67, 93, 111, 111, 62, 62, 111, 111, 36, 36, 
    49, 155, 49, 62, 49, 49, 62, 62, 10, 36, 36, 36, 124, 124, 124, 36, 124
]

find_next_number(sequence)
