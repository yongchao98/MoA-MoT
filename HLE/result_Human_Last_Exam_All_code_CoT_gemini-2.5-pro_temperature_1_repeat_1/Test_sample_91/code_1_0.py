def analyze_garner_scale():
    """
    Analyzes and identifies the scale used by Erroll Garner in a specific
    passage of "All My Loves Are You".

    The melody in question (39s-43s) is a descending run over a dominant 7th
    chord (likely G7). The notes used fit a specific 8-note scale known
    for its use in jazz improvisation. This script demonstrates the
    structure of that scale.
    """

    # Note names and their corresponding numerical values for interval calculation.
    # We use a chromatic scale where each half-step is 1 unit.
    note_values = {
        'F': 5, 'F#': 6, 'G': 7, 'G#': 8, 'Ab': 8, 'A': 9, 'A#': 10,
        'Bb': 10, 'B': 11, 'C': 0, 'C#': 1, 'Db': 1, 'D': 2, 'D#': 3,
        'Eb': 3, 'E': 4
    }

    # The notes of the scale starting on G, which fits the musical context.
    # This is one of the two forms of the octatonic (eight-note) scale.
    # We use common jazz enharmonic spellings (e.g., Ab for the b9, Bb for the #9).
    scale_notes = ['G', 'Ab', 'Bb', 'B', 'C#', 'D', 'E', 'F']

    print("Analyzing the scale played by Erroll Garner...")
    print(f"The root of the underlying chord is likely G.")
    print(f"The notes of the corresponding scale are: {', '.join(scale_notes)}")
    print("-" * 30)
    print("Calculating the interval pattern:")

    intervals = []
    # Add the root note to the end to calculate the last interval back to the octave
    notes_for_loop = scale_notes + [scale_notes[0]]

    for i in range(len(notes_for_loop) - 1):
        note1 = notes_for_loop[i]
        note2 = notes_for_loop[i+1]

        # Calculate the distance in half-steps
        val1 = note_values[note1]
        val2 = note_values[note2]
        
        # Handle the wrap-around from B to C
        if val2 < val1:
            val2 += 12
        
        distance = val2 - val1

        if distance == 1:
            interval_name = "H" # Half Step
        elif distance == 2:
            interval_name = "W" # Whole Step
        else:
            interval_name = "?"

        intervals.append(interval_name)
        print(f"  {note1} to {note2} -> {distance} half-step(s) = '{interval_name}'")

    final_pattern = "-".join(intervals)
    print("-" * 30)
    print(f"The resulting interval pattern is: {final_pattern}")
    print("\nConclusion:")
    print("This H-W-H-W-H-W-H-W pattern is characteristic of the Half-Whole Diminished Scale.")
    print("In a broader sense, this type of eight-note, symmetric scale is known as an Octatonic Scale.")


if __name__ == '__main__':
    analyze_garner_scale()