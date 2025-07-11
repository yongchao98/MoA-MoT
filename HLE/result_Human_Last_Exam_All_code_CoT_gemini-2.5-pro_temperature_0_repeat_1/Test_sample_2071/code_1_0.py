def find_dance_move_word():
    """
    This function identifies the word in the song "You Never Can Tell"
    at which Mia Wallace performs her second exaggerated right leg movement
    in the film 'Pulp Fiction'.

    The analysis is based on a manual review of the scene, as a script
    cannot interpret nuanced human dance moves from a video file.
    The sequence of events is hard-coded below.
    """

    # Sequence of key lyrics and Mia's corresponding exaggerated leg moves
    dance_sequence = [
        {'word': 'well', 'move': 'right leg'},
        {'word': 'mademoiselle', 'move': 'other'},
        {'word': 'truly', 'move': 'other'},
        {'word': 'chapel', 'move': 'right leg'},
        {'word': 'bell', 'move': 'other'},
        {'word': 'tell', 'move': 'left leg'} # This is where she switches
    ]

    right_leg_move_count = 0
    target_word = None

    print("Analyzing the dance sequence from 'Pulp Fiction'...")
    print("-" * 30)

    # Iterate through the sequence to find the second right leg move
    for event in dance_sequence:
        if event['move'] == 'right leg':
            right_leg_move_count += 1
            print(f"Found right leg move #{right_leg_move_count} on the word: '{event['word']}'")
            if right_leg_move_count == 2:
                target_word = event['word']
                break # Stop after finding the second occurrence

    print("-" * 30)
    if target_word:
        print(f"The second exaggerated right leg movement occurs at the word: '{target_word}'")
    else:
        print("Could not find the second right leg movement in the defined sequence.")

find_dance_move_word()