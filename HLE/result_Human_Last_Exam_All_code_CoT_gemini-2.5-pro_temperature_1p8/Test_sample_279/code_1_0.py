def get_piano_fingering():
    """
    Determines and prints the most comfortable fingering for a wide piano chord.
    
    The chord is too wide to be played at once, so it's treated as a roll/arpeggio.
    The fingering uses a standard thumb-cross-under technique for comfort and fluidity.
    """
    
    # The notes of the chord
    notes = ["Bb3", "Db4", "F4", "Bb4", "Db5"]
    
    # The most comfortable fingering for playing this as an arpeggio
    fingering = ["1", "2", "3", "1", "2"]
    
    # Combine the notes and fingerings into the desired output format
    # Example: "Bb3 1", "Db4 2", etc.
    formatted_parts = []
    for i in range(len(notes)):
        part = f"{notes[i]} {fingering[i]}"
        formatted_parts.append(part)
        
    # Join the parts with a comma and space
    final_output = ", ".join(formatted_parts)
    
    print(final_output)

get_piano_fingering()
<<<Bb3 1, Db4 2, F4 3, Bb4 1, Db5 2>>>