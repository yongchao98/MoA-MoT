def solve_sequence():
    """
    This function generates the sequence based on the identified pattern
    and prints the next 4 elements.
    """
    # The initial sequence provided by the user.
    initial_sequence = [3, 2, 1, 2, 3, 3, 3, 2, 2]

    # The pattern for the values of the runs is a repeating cycle of (3, 2, 1, 2).
    run_values_pattern = [3, 2, 1, 2]
    
    # The pattern for the lengths of the runs cycles through (1,1,1,1) and then (3,3,3,3).
    run_lengths_pattern = [1, 1, 1, 1, 3, 3, 3, 3]

    # We generate the full sequence programmatically to find the next elements.
    generated_sequence = []
    run_index = 0
    # Generate enough elements to find the next 4.
    while len(generated_sequence) < len(initial_sequence) + 4:
        # Get the value and length for the current run from the patterns.
        value = run_values_pattern[run_index % len(run_values_pattern)]
        length = run_lengths_pattern[run_index % len(run_lengths_pattern)]
        
        # Add the run to the sequence.
        for _ in range(length):
            generated_sequence.append(value)
        
        run_index += 1

    # The full sequence including the next 4 elements.
    final_sequence = generated_sequence[:len(initial_sequence) + 4]

    # The instruction is to output each number in the final equation.
    # We will print the full sequence including the newly found elements.
    print("The final sequence is:")
    print(' '.join(map(str, final_sequence)))

solve_sequence()