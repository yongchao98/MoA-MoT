def solve():
    """
    Determines and prints the most comfortable right-hand piano fingering for the given chord.
    The chord is Bb3, Db4, F4, Bb4, Db5.
    Given the wide span, it must be played as a roll.
    """
    
    # The notes of the chord in order from lowest to highest.
    notes = ["Bb3", "Db4", "F4", "Bb4", "Db5"]
    
    # The most comfortable fingering for a right-hand roll.
    # 1 (thumb) on the lowest note.
    # 2 (index), 3 (middle), 4 (ring), 5 (pinky) on the upper notes.
    fingers = [1, 2, 3, 4, 5]
    
    # Build the output string in the format "Note Finger, Note Finger, ..."
    result_parts = []
    for i in range(len(notes)):
        note = notes[i]
        finger = fingers[i]
        # The problem statement requires printing each number in the final equation.
        # This format "Note Finger" satisfies that.
        result_parts.append(f"{note} {finger}")
        
    final_answer = ", ".join(result_parts)
    
    print(final_answer)

solve()