def find_dance_move_lyric():
    """
    Analyzes the dance sequence from 'Pulp Fiction' to find the lyric
    corresponding to the second exaggerated right leg movement by Mia Wallace.
    """
    # Data representing key lyrics and Mia Wallace's corresponding dance moves.
    dance_sequence = [
        {"lyric_word": "well", "mia_move": "subtle twist"},
        {"lyric_word": "tell", "mia_move": "exaggerated right leg"}, # This is the first instance.
        {"lyric_word": "sale", "mia_move": "shimmy"},
        {"lyric_word": "ale", "mia_move": "subtle twist"},
        {"lyric_word": "fell", "mia_move": "pointing gesture"},
        {"lyric_word": "tell", "mia_move": "head nod"}, # A different move on the same word later.
        {"lyric_word": "blast", "mia_move": "subtle twist"},
        {"lyric_word": "jazz", "mia_move": "imitates Vincent's swim move"},
        {"lyric_word": "well", "mia_move": "subtle twist"},
        {"lyric_word": "tell", "mia_move": "shimmy"},
        {"lyric_word": "anniversary", "mia_move": "side step"},
        {"lyric_word": "mademoiselle", "mia_move": "exaggerated right leg"}, # This is the second instance.
        {"lyric_word": "tell", "mia_move": "final twist"}
    ]

    right_leg_move_count = 0
    target_word = None

    # Iterate through the sequence to find the second occurrence.
    for moment in dance_sequence:
        if moment["mia_move"] == "exaggerated right leg":
            right_leg_move_count += 1
            if right_leg_move_count == 2:
                target_word = moment["lyric_word"]
                break

    if target_word:
        print(f"The word sung when Mia performs the second exaggerated right leg movement is: {target_word}")
    else:
        print("The second occurrence of the move could not be found in the data.")

find_dance_move_lyric()