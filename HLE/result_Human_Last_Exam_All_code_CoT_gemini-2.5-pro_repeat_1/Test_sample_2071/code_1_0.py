def find_dance_move_lyric():
    """
    Analyzes the dance sequence from 'Pulp Fiction' to find the lyric
    corresponding to a specific dance move.
    """
    # This list represents the chronological sequence of Mia Wallace's key dance
    # moves synchronized with the lyrics of "You Never Can Tell".
    # Each entry is a tuple: (sung_word, mia_primary_move).
    scene_analysis = [
        ("well", "General twisting"),
        ("truly", "Starts first right leg exaggeration"),
        ("mademoiselle", "Continues right leg exaggeration"),
        ("chapel", "Switches to left leg exaggeration"),
        ("bell", "Continues left leg exaggeration"),
        ("tell", "Switches to second right leg exaggeration")
    ]

    right_leg_move_count = 0
    target_word = None

    # We iterate through the sequence to find the second instance
    # of a move involving the right leg.
    for word, move in scene_analysis:
        if "right leg exaggeration" in move:
            right_leg_move_count += 1
            if right_leg_move_count == 2:
                target_word = word
                break

    if target_word:
        print(f"In the song 'You Never Can Tell', the dance scene analysis shows that:")
        print(f"1. The first right leg exaggeration begins at the word: 'truly'")
        print(f"2. The switch to the left leg occurs at the word: 'chapel'")
        print(f"3. The second right leg exaggeration begins at the word: '{target_word}'")
    else:
        print("Could not determine the specific word from the analysis.")

find_dance_move_lyric()