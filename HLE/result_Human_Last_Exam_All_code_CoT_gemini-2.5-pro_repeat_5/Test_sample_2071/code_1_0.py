def find_dance_move_lyric():
    """
    This function analyzes the dance sequence from 'Pulp Fiction' to find the
    lyric corresponding to the second exaggerated right leg move.
    """
    # This data represents the timeline of key dance moves and corresponding lyrics.
    # Each tuple contains: (Lyric Word, Description of Mia's primary move)
    dance_sequence = [
        ("well", "First exaggerated right leg kick out."),
        ("mademoiselle", "Shift to exaggerated left leg movement."),
        ("bell", "Continued focus on left leg movement."),
        ("tell", "First time for the 'hands over eyes' move."),
        ("sale", "Side-to-side steps."),
        ("ale", "Side-to-side steps."),
        ("well", "Second exaggerated right leg kick out."),
        ("tell", "Second time for the 'hands over eyes' move."),
        ("blast", "Twisting moves."),
        ("jazz", "Twisting moves."),
        ("fell", "Sliding move across the floor."),
        ("tell", "Third time for the 'hands over eyes' move.")
    ]

    target_move = "exaggerated right leg kick"
    target_occurrence = 2
    current_occurrence = 0
    found_lyric = None

    print(f"Searching for the lyric at occurrence number {target_occurrence} of the move: '{target_move}'")
    print("--------------------------------------------------")

    for i, (lyric, move) in enumerate(dance_sequence):
        # Check if the move description contains our target phrase
        if target_move in move:
            current_occurrence += 1
            print(f"Found occurrence {current_occurrence} at lyric: '{lyric}'")
            if current_occurrence == target_occurrence:
                found_lyric = lyric
                break
    
    print("--------------------------------------------------")
    if found_lyric:
        # Final "equation" style output
        # Lyric at (Move Occurrence 2) = Result
        print(f"The lyric sung during the second exaggerated right leg movement is:")
        print(f"Lyric at (Move Occurrence {target_occurrence}) = {found_lyric}")
    else:
        print("Could not find the specified occurrence of the move.")

find_dance_move_lyric()