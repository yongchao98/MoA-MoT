def get_piano_fingering():
    """
    This function determines and prints the most comfortable right-hand fingering
    for a wide arpeggiated chord.
    """
    # The notes of the chord in order
    notes = ["Bb3", "Db4", "F4", "Bb4", "Db5"]
    
    # The most comfortable fingering for playing these notes as a right-hand arpeggio.
    # 1: Thumb, 2: Index, 3: Middle, 4: Ring, 5: Pinky
    # The fingering is 1-2-3, then a thumb-tuck to 1, followed by 2.
    fingering = [1, 2, 3, 1, 2]
    
    # Build the output string in the format "Note Finger, Note Finger, ..."
    output_parts = []
    for i in range(len(notes)):
        part = f"{notes[i]} {fingering[i]}"
        output_parts.append(part)
        
    final_output = ", ".join(output_parts)
    
    print(final_output)

get_piano_fingering()