def get_piano_fingering():
    """
    Determines and prints the most comfortable right-hand fingering for the specified chord.
    """
    # The notes of the chord in ascending order.
    notes = ["Bb3", "Db4", "F4", "Bb4", "Db5"]
    
    # The standard and most comfortable fingering for this chord shape.
    # 1: Thumb, 2: Index, 3: Middle, 4: Ring, 5: Pinky
    fingering = [1, 2, 3, 4, 5]
    
    # Prepare the output strings for each note-finger pair.
    output_parts = []
    for i in range(len(notes)):
        note = notes[i]
        finger = fingering[i]
        # The final output needs to include each number in the equation.
        # So we format it as "Note Finger".
        output_parts.append(f"{note} {finger}")
        
    # Join all parts into a single string, separated by ", ".
    final_output_string = ", ".join(output_parts)
    
    print(final_output_string)

get_piano_fingering()
<<<Bb3 1, Db4 2, F4 3, Bb4 4, Db5 5>>>