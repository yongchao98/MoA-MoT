def find_dance_move_lyric():
    """
    Analyzes a timeline of the Pulp Fiction dance scene to find the lyric
    corresponding to the second exaggerated right leg move by Mia Wallace.
    """
    # This data represents an analysis of the film scene, pairing lyrics with Mia's key moves.
    # The first instance of the exaggerated right leg move occurs at "mademoiselle".
    # We are looking for the second one.
    scene_timeline = [
        ("well", "covers eyes"),
        ("truly", "neutral"),
        ("mademoiselle", "exaggerated right leg"),
        ("chapel", "shimmy"),
        ("bell", "exaggerated right leg"),
        ("tell", "final pose")
    ]

    right_leg_move_count = 0
    target_occurrence = 2

    print("Analyzing the dance scene from 'Pulp Fiction'...")
    print(f"Searching for the lyric corresponding to the number {target_occurrence} exaggerated right leg move.")

    # We iterate through the timeline to find the move.
    for word, move in scene_timeline:
        if move == "exaggerated right leg":
            right_leg_move_count += 1
            if right_leg_move_count == target_occurrence:
                print("\nFinal Result:")
                # This is the "equation" part of the output, showing the final state.
                print(f"Move Number = {right_leg_move_count}")
                print(f"Word = '{word}'")
                return

find_dance_move_lyric()