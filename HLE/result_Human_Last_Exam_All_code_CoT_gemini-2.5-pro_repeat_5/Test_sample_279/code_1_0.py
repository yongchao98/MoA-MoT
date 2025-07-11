def get_piano_fingering():
    """
    Determines and prints the most comfortable right-hand fingering for a given chord.

    The chord consists of the notes: Bb3, Db4, F4, Bb4, Db5.
    The fingering is determined based on ergonomic principles for piano playing.
    """
    # The notes of the chord in order from lowest to highest.
    notes = ['Bb3', 'Db4', 'F4', 'Bb4', 'Db5']

    # The most comfortable and logical fingering for a wide chord roll on the right hand.
    # 1 = Thumb, 2 = Index, 3 = Middle, 4 = Ring, 5 = Pinky.
    # This sequential fingering (1-2-3-4-5) is the most standard and comfortable
    # for playing wide, arpeggiated or rolled chords like this one.
    fingers = [1, 2, 3, 4, 5]

    # Create a list to hold the formatted "Note Finger" strings.
    formatted_fingering = []
    for i in range(len(notes)):
        note = notes[i]
        finger = fingers[i]
        # The final equation includes each number for each note
        formatted_fingering.append(f"{note} {finger}")

    # Join the list elements with a comma and space for the final output.
    final_answer = ", ".join(formatted_fingering)
    
    print(final_answer)
    print(f"<<<{final_answer}>>>")

get_piano_fingering()