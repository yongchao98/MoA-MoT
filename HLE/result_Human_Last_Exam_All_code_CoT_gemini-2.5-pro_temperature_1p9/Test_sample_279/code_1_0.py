def solve_fingering():
    """
    Determines and prints the most comfortable right-hand piano fingering
    for the specified broken chord.
    """
    notes = ["Bb3", "Db4", "F4", "Bb4", "Db5"]
    fingering = [1, 2, 3, 1, 2]
    
    # Combine the notes and fingering into the desired format
    result_parts = []
    for i in range(len(notes)):
        result_parts.append(f"{notes[i]} {fingering[i]}")
        
    final_string = ", ".join(result_parts)
    
    print(final_string)

solve_fingering()