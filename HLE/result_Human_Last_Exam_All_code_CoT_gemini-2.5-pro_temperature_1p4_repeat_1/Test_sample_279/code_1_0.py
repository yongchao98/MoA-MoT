def get_piano_fingering():
    """
    Determines and prints the most comfortable right-hand fingering for a given piano chord.
    """
    notes = ["Bb3", "Db4", "F4", "Bb4", "Db5"]
    # The most comfortable and standard fingering for this chord is 1-2-3-4-5.
    fingers = [1, 2, 3, 4, 5]

    # Combine notes and fingerings into the desired format.
    # The format is: Note1 Finger1, Note2 Finger2, ...
    result_parts = []
    for i in range(len(notes)):
        note = notes[i]
        finger = fingers[i]
        # We need to print each number in the final equation.
        # Here we format each part like "Bb3 1"
        result_parts.append(f"{note} {finger}")

    # Join all the parts with a comma and a space.
    final_string = ", ".join(result_parts)

    print(final_string)

get_piano_fingering()
print("<<<Bb3 1, Db4 2, F4 3, Bb4 4, Db5 5>>>")