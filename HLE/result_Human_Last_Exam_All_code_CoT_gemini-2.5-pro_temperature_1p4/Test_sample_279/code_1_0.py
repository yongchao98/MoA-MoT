def get_piano_fingering():
    """
    Determines and prints the most comfortable right-hand fingering for a wide-spread chord.

    The logic is based on standard piano ergonomic principles for wide chords:
    1. The chord is Bb3, Db4, F4, Bb4, Db5. The span is a minor 12th, which is very wide.
    2. The most logical fingering for a hand large enough to play this is a sequential one.
    3. The thumb (1) takes the lowest note, and the pinky (5) takes the highest.
    4. The inner fingers (2, 3, 4) take the inner notes in order to avoid awkward crossovers.
    """
    
    # The notes of the chord in ascending order
    notes = ['Bb3', 'Db4', 'F4', 'Bb4', 'Db5']
    
    # The corresponding finger for each note (1=thumb, 5=pinky)
    # This sequential fingering is the most direct and standard, though it requires a large hand.
    fingerings = [1, 2, 3, 4, 5]

    # Build the list of "note finger" strings
    output_parts = []
    for i in range(len(notes)):
        note = notes[i]
        finger = fingerings[i]
        # Each part of the final output string is formatted here
        # This uses the numbers from the notes (e.g., 3, 4, 4, 4, 5) and the fingerings (1, 2, 3, 4, 5)
        output_parts.append(f"{note} {finger}")

    # Join the parts with a comma and a space
    final_answer = ", ".join(output_parts)
    
    # Print the final result
    print(final_answer)

get_piano_fingering()