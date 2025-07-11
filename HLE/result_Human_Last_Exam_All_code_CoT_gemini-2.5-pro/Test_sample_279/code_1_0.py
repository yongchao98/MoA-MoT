def get_piano_fingering():
    """
    This function determines and prints the most comfortable right-hand fingering
    for the chord Bb3, Db4, F4, Bb4, Db5.

    The chord is very wide and must be rolled. The fingering 1-2-3-4-5
    is the most ergonomic for this motion.
    """
    notes = ["Bb3", "Db4", "F4", "Bb4", "Db5"]
    fingers = [1, 2, 3, 4, 5]

    # Build the final output string by pairing each note with its finger number.
    output_parts = []
    for i in range(len(notes)):
        note = notes[i]
        finger = fingers[i]
        output_parts.append(f"{note} {finger}")

    final_answer = ", ".join(output_parts)
    print(final_answer)

# Execute the function to display the result.
get_piano_fingering()