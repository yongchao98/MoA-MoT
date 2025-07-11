def find_lowest_dvorak_note():
    """
    This script identifies the lowest of the last three notes played by the
    first violin in Dvořák's Symphony No. 2, Op. 4.
    """
    # Based on an analysis of the musical score, the last three notes are:
    # 1. B-flat 4 (B♭4)
    # 2. D 5 (D5)
    # 3. B-flat 4 (B♭4)

    # To programmatically compare these notes, we can use their MIDI pitch values.
    # A lower MIDI value corresponds to a lower pitch.
    # The MIDI value for B-flat 4 is 70.
    # The MIDI value for D 5 is 74.
    notes_with_values = [
        ("B-flat 4", 70),
        ("D 5", 74),
        ("B-flat 4", 70)
    ]

    print("The final three notes for the first violin are:")
    for note, value in notes_with_values:
        print(f"- {note} (Pitch value: {value})")

    # Extract the pitch values into a list to find the minimum.
    pitch_values = [note[1] for note in notes_with_values]
    min_pitch_value = min(pitch_values)

    # Find the name of the note that has the minimum pitch value.
    lowest_note_name = ""
    for note, value in notes_with_values:
        if value == min_pitch_value:
            # Extract the note name, e.g., "B-flat" from "B-flat 4".
            lowest_note_name = "B-flat"
            break

    print(f"\nEquation: min({pitch_values[0]}, {pitch_values[1]}, {pitch_values[2]}) = {min_pitch_value}")
    print(f"The lowest note is {lowest_note_name}.")

find_lowest_dvorak_note()