def get_piano_fingering():
    """
    Determines and prints the most comfortable fingering for a given piano chord.
    The final output must include each note of the chord in the final equation.
    """
    notes = ['Bb3', 'Db4', 'F4', 'Bb4', 'Db5']
    fingering = [1, 2, 3, 4, 5]

    # Combine the notes and fingering into the desired format
    result_parts = []
    for i in range(len(notes)):
        # Construct the string for each note, e.g., "Bb3 1"
        part = f"{notes[i]} {fingering[i]}"
        result_parts.append(part)

    # Join all parts with a comma and a space
    final_output = ", ".join(result_parts)

    print(final_output)

get_piano_fingering()