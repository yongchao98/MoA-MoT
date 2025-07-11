def get_piano_fingering():
    """
    Determines and prints the most comfortable right-hand piano fingering
    for the specified chord: Bb3, Db4, F4, Bb4, Db5.
    """
    notes = ["Bb3", "Db4", "F4", "Bb4", "Db5"]
    fingering = [1, 2, 3, 4, 5]

    # Create the formatted string for the final answer
    # e.g., "Bb3 1, Db4 2, F4 3, Bb4 4, Db5 5"
    output_parts = []
    for i in range(len(notes)):
        part = f"{notes[i]} {fingering[i]}"
        output_parts.append(part)
    
    final_answer = ", ".join(output_parts)
    
    print("The most comfortable fingering for the chord is:")
    print(final_answer)
    print("\n<<<" + final_answer + ">>>")

get_piano_fingering()