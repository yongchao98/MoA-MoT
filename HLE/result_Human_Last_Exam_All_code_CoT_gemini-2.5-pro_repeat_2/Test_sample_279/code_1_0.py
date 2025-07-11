def get_piano_fingering():
    """
    Determines and prints the most comfortable right-hand piano fingering
    for the specified chord.
    """

    # The notes and their corresponding finger numbers based on pianistic principles.
    # 1: Thumb, 2: Index, 3: Middle, 4: Ring, 5: Pinky
    chord_fingering = {
        "Bb3": 1,
        "Db4": 2,
        "F4": 3,
        "Bb4": 4,
        "Db5": 5
    }

    # Format the fingering into the desired string format.
    # Example: "Note1 Finger1, Note2 Finger2, ..."
    output_string = ", ".join([f"{note} {finger}" for note, finger in chord_fingering.items()])

    print(output_string)

get_piano_fingering()