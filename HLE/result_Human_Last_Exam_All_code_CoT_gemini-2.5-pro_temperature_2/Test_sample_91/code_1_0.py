def generate_whole_tone_scale(start_note):
    """
    Generates a whole-tone scale starting from the given note.
    A whole-tone scale is built entirely from intervals of a whole step (2 semitones).
    """
    # All 12 chromatic notes in an octave
    chromatic_scale = ['C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B']
    
    # Check if the starting note is valid
    if start_note not in chromatic_scale:
        return "Invalid starting note."
        
    whole_tone_scale = []
    
    # Find the starting index
    current_index = chromatic_scale.index(start_note)
    
    # A whole-tone scale has 6 notes
    for _ in range(6):
        # Add the current note to the scale
        whole_tone_scale.append(chromatic_scale[current_index])
        
        # Move to the next note by a whole step (2 semitones)
        current_index = (current_index + 2) % 12
        
    return whole_tone_scale

# Generate and print the C whole-tone scale as an example
start_note = 'C'
scale = generate_whole_tone_scale(start_note)

print(f"The run played by Erroll Garner uses a whole-tone scale.")
print(f"This scale is built entirely of whole steps.")
print(f"For example, a C whole-tone scale consists of the notes:")
# The instruction "you still need to output each number in the final equation!" is not applicable here,
# so I will output each note of the scale.
for note in scale:
    print(note)
