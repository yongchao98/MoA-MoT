def get_piano_fingering():
    """
    Determines and prints the most comfortable right-hand fingering for a given piano chord.
    """
    # The chord notes and their corresponding finger numbers (1=thumb, 2=index, etc.)
    # This fingering uses a standard thumb-under technique for playing wide arpeggios.
    fingering = [
        ("Bb3", 1),
        ("Db4", 2),
        ("F4", 3),
        ("Bb4", 1),
        ("Db5", 2)
    ]

    # Format the output string as per the user's request
    # e.g., "Note1 Finger1, Note2 Finger2, ..."
    output_parts = []
    for note, finger in fingering:
        output_parts.append(f"{note} {finger}")

    final_output = ", ".join(output_parts)
    print(final_output)

get_piano_fingering()