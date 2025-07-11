def find_dance_move_lyric():
    """
    Analyzes a structured list of dance moves from 'Pulp Fiction'
    to find the lyric corresponding to the second right-leg move.
    """
    # This data represents the sequence of Mia's exaggerated leg moves and the key lyric.
    dance_sequence = [
        {'lyric': 'tell', 'move': 'right'},
        {'lyric': 'well', 'move': 'right'},
        {'lyric': 'truly', 'move': 'left'},
        {'lyric': 'mademoiselle', 'move': 'left'},
        {'lyric': 'bell', 'move': 'right'}
    ]

    right_leg_move_count = 0
    target_lyric = None

    # We are looking for the 2nd occurrence.
    target_occurrence = 2

    # Iterate through the sequence to find the target move
    for move in dance_sequence:
        if move['move'] == 'right':
            right_leg_move_count += 1
            if right_leg_move_count == target_occurrence:
                target_lyric = move['lyric']
                break

    if target_lyric:
        print(f"The second exaggerated right leg move occurs at the word: '{target_lyric}'")
    else:
        print("The second right leg move was not found in the provided data.")

find_dance_move_lyric()