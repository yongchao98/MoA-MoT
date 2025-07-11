def generate_whole_tone_scale(start_note):
    """
    Generates a 6-note whole-tone scale starting from a given note.
    A whole-tone scale is built entirely from whole steps (2 semitones).
    """
    chromatic_scale = ['C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B']
    
    # Find the starting index
    try:
        current_index = chromatic_scale.index(start_note.upper())
    except ValueError:
        return f"Error: Note '{start_note}' not found in the chromatic scale."

    whole_tone_scale = []
    # A whole-tone scale has 6 notes
    for _ in range(6):
        # Add the current note to the scale
        whole_tone_scale.append(chromatic_scale[current_index])
        # Move to the next note by a whole step (2 semitones)
        current_index = (current_index + 2) % 12
        
    return whole_tone_scale

# Analysis of the melody between 0:39 and 0:43 of Erroll Garner's
# "All My Loves Are You" reveals the use of a whole-tone scale.
# This script will generate the notes of that scale.
scale_type = "Whole-Tone Scale"
starting_note_of_run = "F"
garner_scale = generate_whole_tone_scale(starting_note_of_run)

print(f"The scale Erroll Garner plays in the right hand is a {scale_type}.")
print("The notes in this scale are:")

# Print the final equation (the notes of the scale)
print(' - '.join(garner_scale))
