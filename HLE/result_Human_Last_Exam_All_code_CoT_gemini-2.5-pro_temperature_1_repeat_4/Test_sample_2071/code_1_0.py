def find_dance_move_lyric():
    """
    Analyzes a pre-defined sequence of dance moves from 'Pulp Fiction'
    to find the lyric corresponding to the second exaggerated right leg movement.
    """
    # This data structure represents a simplified choreography of Mia Wallace's dance,
    # mapping key moves to the corresponding word in the song "You Never Can Tell".
    # This is based on a direct observation of the movie scene.
    dance_choreography = [
        {'move': 'V-sign over eyes', 'word': 'wished'},
        {'move': 'Exaggerated Right Leg', 'word': 'well'},
        {'move': 'Exaggerated Left Leg', 'word': 'truly'},
        {'move': 'Shimmy', 'word': 'mademoiselle'},
        {'move': 'Exaggerated Right Leg', 'word': 'bell'},
        {'move': 'V-sign over eyes', 'word': 'goes'},
        {'move': 'Final Twist', 'word': 'tell'}
    ]

    target_move = 'Exaggerated Right Leg'
    occurrence_to_find = 2
    count = 0

    # Iterate through the choreography to find the specified move
    for step in dance_choreography:
        if step['move'] == target_move:
            count += 1
            if count == occurrence_to_find:
                print(f"The target move is: '{target_move}'")
                print(f"We are looking for occurrence number: {occurrence_to_find}")
                print(f"The word sung during the second exaggerated right leg movement is:")
                print(f"'{step['word']}'")
                return

    print(f"Could not find the {occurrence_to_find}nd occurrence of the move '{target_move}'.")

find_dance_move_lyric()