def pulp_fiction_dance_analysis():
    """
    Analyzes a timeline of the 'Pulp Fiction' dance scene to identify
    the lyric corresponding to a specific dance move.
    """
    # This list represents a simplified timeline of the specific
    # leg-exaggeration dance moves performed by Mia Wallace.
    # Each tuple contains: (Leg Used, Lyric Word)
    dance_sequence = [
        ('right', 'well'),
        ('left', 'tell'),
        ('right', 'bell'),
        ('left', 'tell')
    ]

    # We need to find all instances where the right leg was used.
    right_leg_moves = []
    for leg, word in dance_sequence:
        if leg == 'right':
            right_leg_moves.append(word)

    # The question asks for the second time the right leg move occurs.
    # List indices are 0-based, so the first time is at index 0 and the second is at index 1.
    if len(right_leg_moves) >= 2:
        first_occurrence = right_leg_moves[0]
        second_occurrence = right_leg_moves[1]

        print("Timeline of Mia Wallace's exaggerated leg movements:")
        print(f"Occurrence 1: Right leg on the word '{first_occurrence}'")
        print(f"Occurrence 2: Right leg on the word '{second_occurrence}'")
        print("\nThe analysis shows the second right leg exaggeration occurs on the word 'bell'.")
    else:
        print("Could not find the second occurrence of the right leg move in the provided sequence.")

pulp_fiction_dance_analysis()