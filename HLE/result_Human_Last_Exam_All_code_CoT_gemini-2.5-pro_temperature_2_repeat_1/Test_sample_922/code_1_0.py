def solve_sequence():
    """
    This function solves for the next number in a sequence of known integer values.
    The sequence is not mathematical but based on a series of discoveries.
    """

    # The sequence of the first six known integer values.
    sequence = [24663, 35005, 119261, 196219, 211770, 227296]

    # As of August 2022, the single known integer that completes this sequence was discovered.
    next_value = 241511

    # The completed sequence.
    completed_sequence = sequence + [next_value]
    
    # Per the instruction to "output each number in the final equation",
    # we will display the progression to the final answer.
    equation_str = " -> ".join(map(str, completed_sequence))

    print("The completed sequence is:")
    print(equation_str)
    
    print("\nThe next number in the sequence is:")
    print(next_value)

solve_sequence()