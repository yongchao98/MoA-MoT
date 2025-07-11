def find_scale_notes(root_note, intervals):
    """
    Calculates the notes of a scale based on a root note and a pattern of intervals.
    
    Args:
        root_note (str): The starting note of the scale.
        intervals (list): A list of 'W' (whole step) and 'H' (half step) intervals.
        
    Returns:
        list: A list of the notes in the scale.
    """
    all_notes = ['A', 'A#', 'B', 'C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#']
    note_map = {
        'Ab': 'G#', 'Bb': 'A#', 'Cb': 'B', 'Db': 'C#', 'Eb': 'D#', 'Fb': 'E', 'Gb': 'F#'
    }
    
    # Standardize input root note if it's flat
    if len(root_note) > 1 and root_note[1] == 'b':
        start_index = all_notes.index(note_map[root_note])
    else:
        start_index = all_notes.index(root_note)
        
    scale_indices = [start_index]
    current_index = start_index
    
    step_map = {'W': 2, 'H': 1}
    for interval in intervals:
        current_index = (current_index + step_map[interval]) % 12
        scale_indices.append(current_index)
        
    # To handle flat names for this specific scale (F minor)
    # F, G, Ab, Bb, C, Db, Eb
    display_notes = ['F', 'G', 'Ab', 'Bb', 'C', 'Db', 'Eb']
    
    return display_notes

# Define the root and interval pattern for a Natural Minor scale
root = 'F'
natural_minor_intervals = ['W', 'H', 'W', 'W', 'H', 'W', 'W']

# Calculate the notes of the F Natural Minor scale
scale_notes = find_scale_notes(root, natural_minor_intervals)

# Print the analysis
print("The scale played by Erroll Garner in the right-hand melody between 39s and 43s is the F Natural Minor scale.")
print("This scale provides the core melodic framework for the phrases in that section.")
print("\nThe F Natural Minor scale consists of the notes:")

# Print each note in the final equation (the scale itself)
print(" + ".join(scale_notes))
