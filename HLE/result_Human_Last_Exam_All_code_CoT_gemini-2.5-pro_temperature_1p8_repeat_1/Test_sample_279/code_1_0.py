def get_piano_fingering():
    """
    Determines and formats the most comfortable piano fingering for a given wide chord.
    """
    # The notes of the chord provided by the user.
    notes = ["Bb3", "Db4", "F4", "Bb4", "Db5"]
    
    # The most logical and comfortable fingering for this rolled chord.
    # 1=Thumb, 2=Index, 3=Middle, 4=Ring, 5=Pinky.
    fingering = [1, 2, 3, 4, 5]
    
    # We will build a list of strings, e.g., ["Bb3 1", "Db4 2", ...]
    formatted_parts = []
    
    # Loop through the notes and fingerings to create the formatted parts.
    # The 'zip' function pairs each note with its corresponding finger number.
    for note, finger in zip(notes, fingering):
        formatted_parts.append(f"{note} {finger}")
        
    # Join the parts together with a comma and a space.
    final_output = ", ".join(formatted_parts)
    
    # Print the final result.
    print(final_output)

get_piano_fingering()
print("<<<Bb3 1, Db4 2, F4 3, Bb4 4, Db5 5>>>")