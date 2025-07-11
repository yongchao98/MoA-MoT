def get_piano_fingering():
    """
    Determines and prints the most comfortable fingering for a wide piano chord.
    The chord is assumed to be rolled, as it's too wide to be played simultaneously.
    """
    notes = ["Bb3", "Db4", "F4", "Bb4", "Db5"]
    # For a rolled chord of this span, the most comfortable fingering is sequential.
    fingers = [1, 2, 3, 4, 5]

    # Create the formatted string
    fingering_parts = []
    for i in range(len(notes)):
        fingering_parts.append(f"{notes[i]} {fingers[i]}")

    final_fingering = ", ".join(fingering_parts)
    print(final_fingering)

get_piano_fingering()