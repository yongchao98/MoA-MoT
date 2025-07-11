def solve_sequence():
    """
    This function solves the sequence puzzle by identifying a structural pattern.
    """
    initial_sequence = [3, 2, 1, 2, 3, 3, 3, 2, 2]

    # The pattern identified is based on grouping the numbers:
    # Group 1: [3, 2, 1] -> A descending sequence.
    # Group 2: [2, 3] -> An ascending sequence.
    # Group 3: [3, 3, 3] -> A block of three '3's.
    # Group 4: [2, 2, ...] -> The start of a block of '2's.
    #
    # Logic: The pattern suggests that after a block of three '3's, there will be a
    # block of three '2's, followed by a block of three '1's.

    # 1. Complete the block of '2's.
    next_element_1 = 2

    # 2. Add the next block, which is three '1's.
    next_element_2 = 1
    next_element_3 = 1
    next_element_4 = 1

    next_four_elements = [next_element_1, next_element_2, next_element_3, next_element_4]

    print("The original sequence is:")
    # Printing each number of the sequence as requested
    print(*initial_sequence)

    print("\nBased on the pattern of descending/ascending groups followed by blocks of numbers,")
    print("the next four elements are determined to be:")
    print(*next_four_elements)

    # The final "equation" is the completed sequence.
    final_sequence = initial_sequence + next_four_elements
    print("\nThe full sequence is:")
    print(*final_sequence)

solve_sequence()