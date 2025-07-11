def get_piano_fingering():
    """
    Determines and prints the most comfortable right-hand fingering for a wide chord.
    The chord must be rolled/arpeggiated.
    """
    notes = ["Bb3", "Db4", "F4", "Bb4", "Db5"]
    # Standard fingering for a wide 5-note arpeggio is sequential.
    fingers = [1, 2, 3, 4, 5]

    # Build the output string in the requested format.
    output_parts = []
    for i in range(len(notes)):
        note = notes[i]
        finger = fingers[i]
        # This part ensures we output each number in the final string as requested.
        output_parts.append(f"{note} {finger}")

    final_fingering = ", ".join(output_parts)
    print(final_fingering)

get_piano_fingering()