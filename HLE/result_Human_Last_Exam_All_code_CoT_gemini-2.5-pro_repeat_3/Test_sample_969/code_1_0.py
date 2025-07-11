def solve_sequence():
    """
    This function solves the sequence puzzle based on the identified pattern.
    """
    # The initial sequence provided by the user
    initial_sequence = [3, 2, 1, 2, 3, 3, 3, 2, 2]

    # The rule is a prefix followed by a generated part.
    prefix_len = 4
    repetition_count = 3
    base_sequence = [3, 2, 1, 2]

    # Create a copy to extend
    full_sequence = list(initial_sequence)
    next_four_elements = []

    # Generate the next 4 elements
    for _ in range(4):
        # Determine the current position in the generated part of the sequence
        num_generated_so_far = len(full_sequence) - prefix_len
        
        # Determine which block of repetitions we are in
        block_index = num_generated_so_far // repetition_count
        
        # Determine which element from the base sequence to use
        base_element_index = block_index % len(base_sequence)
        
        # The next number in the sequence
        next_element = base_sequence[base_element_index]
        
        # Add the new element to our lists
        full_sequence.append(next_element)
        next_four_elements.append(next_element)

    # Print the final sequence as requested
    print("The completed sequence is:")
    print(*full_sequence)
    
    # You can also uncomment the line below to see just the next four elements
    # print("\nThe next 4 elements are:", *next_four_elements)

solve_sequence()