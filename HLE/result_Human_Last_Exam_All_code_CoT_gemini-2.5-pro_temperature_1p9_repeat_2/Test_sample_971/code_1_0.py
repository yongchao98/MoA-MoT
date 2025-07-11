def calculate_turns():
    """
    Calculates the total number of turns performed by Kitri.

    The famous turning sequence in Don Quixote consists of 32 fouettés.
    We can model this as 4 sets of 8 turns.
    """
    musical_phrases = 4
    turns_per_phrase = 8
    total_turns = musical_phrases * turns_per_phrase

    print(f"Based on the structure of the music, with {turns_per_phrase} turns in each of the {musical_phrases} phrases, the total is calculated as:")
    print(f"{musical_phrases} * {turns_per_phrase} = {total_turns}")

calculate_turns()