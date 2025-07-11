def get_note_value(note_name):
    """
    Converts a note in scientific pitch notation (e.g., 'F4', 'C6')
    into a numerical value for easy comparison. A higher value means a higher pitch.
    """
    # A map for the base value of each pitch class.
    pitch_map = {'C': 0, 'D': 2, 'E': 4, 'F': 5, 'G': 7, 'A': 9, 'B': 11}
    
    # Extract the pitch letter and the octave number from the string.
    pitch = note_name[0].upper()
    octave = int(note_name[1:])
    
    # Calculate a comparable value. The formula ensures that higher octaves
    # and higher pitches within an octave result in a larger number.
    value = pitch_map[pitch] + octave * 12
    return value

# Based on the score, the three final notes for the first violin are F4, F5, and C6.
notes = ["F4", "F5", "C6"]

# Use the 'min' function with our custom key to find the note with the lowest pitch value.
lowest_note = min(notes, key=get_note_value)

# As requested, we will output each note as part of a final equation.
# We create a string that represents the operation of finding the lowest note.
equation_str = "Lowest(" + ", ".join(notes) + ") = " + lowest_note

print("The first violin part ends on a chord of three notes.")
print("To find the lowest of these notes, we can solve the following equation:")
print(equation_str)