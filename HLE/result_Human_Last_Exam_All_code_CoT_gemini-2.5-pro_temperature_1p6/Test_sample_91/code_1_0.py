import sys

def generate_scale():
    """
    This function identifies and generates the notes of the scale played by Erroll Garner.
    """
    print("In the song 'All My Loves Are You' (from the album 'Afternoon Of An Elf'), Erroll Garner plays a fast, descending run between 0:39-0:43.")
    print("The scale used is a G Whole-Tone Scale, played over a G7 chord.\n")

    # A reference list of all 12 chromatic notes. We'll use sharps for consistency.
    # Note that D# is the same as Eb, C# is the same as Db, etc.
    ALL_NOTES = ['G', 'G#', 'A', 'A#', 'B', 'C', 'C#', 'D', 'D#', 'E', 'F', 'F#']
    root_note = 'G'
    
    # The whole-tone scale is built with a pattern of 6 whole steps (W).
    # A whole step is an interval of 2 semitones.
    interval_pattern = [2, 2, 2, 2, 2, 2] # in semitones
    
    scale_notes = []
    
    # Find the starting index of the root note
    try:
        current_index = ALL_NOTES.index(root_note)
    except ValueError:
        print(f"Error: Root note '{root_note}' is not valid.", file=sys.stderr)
        return

    # Add the root note to our scale
    scale_notes.append(ALL_NOTES[current_index])

    # Generate the rest of the scale by applying the intervals
    for interval in interval_pattern[:-1]:
        current_index = (current_index + interval) % 12
        scale_notes.append(ALL_NOTES[current_index])

    print("The G Whole-Tone scale contains the following notes (written ascending):")
    print(", ".join(scale_notes))
    print("\nThe 'equation' to build this scale is shown below, where 'W' stands for a Whole Step:")

    # Fulfilling the request to output each note in a final equation format
    print(f"{scale_notes[0]} + W -> {scale_notes[1]} + W -> {scale_notes[2]} + W -> {scale_notes[3]} + W -> {scale_notes[4]} + W -> {scale_notes[5]}")

# Execute the function to print the analysis and the scale
generate_scale()