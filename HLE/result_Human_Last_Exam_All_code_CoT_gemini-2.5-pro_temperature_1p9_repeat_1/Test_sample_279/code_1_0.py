def find_comfortable_fingering():
    """
    Determines and prints a comfortable piano fingering for a wide chord.
    The chord is assumed to be rolled/arpeggiated due to its span.
    """
    notes = ["Bb3", "Db4", "F4", "Bb4", "Db5"]
    # This fingering uses a common arpeggio technique for minor chords (1-2-4)
    # followed by a standard thumb-under cross to handle the wide span.
    # 1. Thumb on the lowest note.
    # 2. 1-2-4 on the first three notes (a relaxed Bb minor arpeggio shape).
    # 3. Thumb crosses under to play the octave Bb.
    # 4. Index finger comfortably plays the final note.
    fingering = [1, 2, 4, 1, 2]

    # Format the output string as requested
    result = []
    for i in range(len(notes)):
        result.append(f"{notes[i]} {fingering[i]}")

    final_output = ", ".join(result)
    print("Because the chord spans more than a tenth, it must be 'rolled' (played as a fast arpeggio).")
    print("By far the most comfortable and standard fingering for this rolled chord is:")
    print(final_output)

find_comfortable_fingering()
<<<Bb3 1, Db4 2, F4 4, Bb4 1, Db5 2>>>