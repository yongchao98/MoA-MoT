def find_lowest_note():
    """
    This script finds the lowest of the last three notes played by the first violin
    in Dvořák's Symphony No. 2, Op. 4.

    The last three unique pitches found in the score for the first violin part are:
    1. B-flat 4 (from a double stop in the penultimate measure)
    2. D 5 (from the same double stop)
    3. B-flat 5 (the final note of the symphony)

    We will convert these notes to MIDI pitch values to determine the lowest one.
    """
    
    # A dictionary mapping the note names to their MIDI pitch values.
    notes = {
        "B-flat 4": 70,
        "D 5": 74,
        "B-flat 5": 82
    }

    print("The last three unique pitches for the first violin are B-flat 4, D 5, and B-flat 5.")
    print("Their corresponding MIDI pitch values are:")
    for note, value in notes.items():
        print(f"- {note}: {value}")
        
    # Create a list of the values for the equation
    pitch_values = list(notes.values())

    # Find the minimum MIDI value
    lowest_pitch_value = min(pitch_values)
    
    # Find the note name that corresponds to the lowest value
    lowest_note_name = ""
    for note, value in notes.items():
        if value == lowest_pitch_value:
            lowest_note_name = note
            break

    # Print the "equation" as requested
    print(f"\nThe equation to find the lowest pitch is: min({pitch_values[0]}, {pitch_values[1]}, {pitch_values[2]})")
    print(f"The result of the equation is: {lowest_pitch_value}")
    print(f"\nThe lowest note is {lowest_note_name}.")

if __name__ == "__main__":
    find_lowest_note()