def solve_sequence():
    """
    This function generates the sequence based on a multi-generational rule
    and prints the next 4 elements.
    """
    # The base pattern of values that repeats.
    base_values = [3, 2, 1, 2]

    # The length of runs for each generation.
    # Gen 1 has run length 1, Gen 2 has run length 3, Gen 3 has run length 4, etc.
    # This pattern (1, 3, 4, 5...) is based on observation of the likely source sequence (OEIS A098155).
    run_lengths = [1, 3, 4, 5, 6]

    # The initial sequence provided by the user.
    initial_sequence = [3, 2, 1, 2, 3, 3, 3, 2, 2]

    # We generate the full sequence until it's long enough to find the next 4 elements.
    full_sequence = []
    generation_index = 0
    # The number of elements we need to find.
    num_to_find = 4
    
    while len(full_sequence) < len(initial_sequence) + num_to_find:
        # Determine the length of runs for the current generation.
        current_run_length = run_lengths[generation_index]
        
        # For each value in the base, create a run of the determined length.
        for value in base_values:
            full_sequence.extend([value] * current_run_length)
        
        generation_index += 1

    # Extract the next 4 elements that follow the initial sequence.
    next_elements = full_sequence[len(initial_sequence):len(initial_sequence) + num_to_find]

    # Print the complete sequence as requested.
    print(" ".join(map(str, initial_sequence)), end=" ")
    print(" ".join(map(str, next_elements)))

solve_sequence()