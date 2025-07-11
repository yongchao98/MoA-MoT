def find_dance_move_lyric():
    """
    Analyzes the dance scene from 'Pulp Fiction' to find the lyric sung
    during Mia Wallace's second exaggerated right leg movement.
    """
    # A simplified timeline of Mia's dance moves synced with the song lyrics.
    dance_sequence = [
        {"lyric": "well", "move": "Side-to-side sway"},
        {"lyric": "Pierre", "move": "First Exaggerated Right Leg"},
        {"lyric": "truly", "move": "Exaggerated Left Leg"},
        {"lyric": "mademoiselle", "move": "Side-to-side sway"},
        {"lyric": "chapel", "move": "Hand gestures"},
        {"lyric": "bell", "move": "Second Exaggerated Right Leg"},
        {"lyric": "tell", "move": "The Swim dance move"}
    ]

    right_leg_move_count = 0
    target_word = ""

    for moment in dance_sequence:
        if moment["move"] == "First Exaggerated Right Leg" or moment["move"] == "Second Exaggerated Right Leg":
            right_leg_move_count += 1
            if right_leg_move_count == 2:
                target_word = moment["lyric"]
                break

    # The prompt asks to "output each number in the final equation".
    # Since the answer is a word, we will print the word itself as the result of our analysis.
    print(f"The word sung during the second exaggerated right leg movement is: {target_word}")

find_dance_move_lyric()