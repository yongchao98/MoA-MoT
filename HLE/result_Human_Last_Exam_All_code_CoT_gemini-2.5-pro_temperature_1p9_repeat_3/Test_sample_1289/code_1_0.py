def solve_music_intervals():
    """
    Calculates and prints the sequence of musical intervals based on the problem description.
    """
    # Step 1: Define the sequence of notes sung by the musician.
    # We represent the notes in terms of semitones from the root note C (C=0).
    # Maqam Zanjaran notes: C(0), D(2), E(4), F(5), G(7), Ab(8), Bb(10), C'(12)
    # The musician ascends C -> C' and then descends C' -> F.
    # The full path of notes sung is:
    # C, D, E, F, G, Ab, Bb, C', Bb, Ab, G, F
    notes_in_semitones = [0, 2, 4, 5, 7, 8, 10, 12, 10, 8, 7, 5]

    # Step 2: Calculate the interval between each consecutive note.
    intervals = []
    for i in range(len(notes_in_semitones) - 1):
        # The interval in semitones is the absolute difference between the note values.
        interval_in_semitones = abs(notes_in_semitones[i+1] - notes_in_semitones[i])
        
        # Convert the interval to the desired unit, where 1.0 is a whole tone (2 semitones).
        interval_in_units = interval_in_semitones / 2.0
        intervals.append(interval_in_units)

    # Step 3: Format the output string as requested.
    # We will format numbers to be integers if they are whole numbers (e.g., 1.0 -> 1)
    # or floats otherwise (e.g., 0.5).
    formatted_intervals = []
    for num in intervals:
        if num == int(num):
            formatted_intervals.append(str(int(num)))
        else:
            formatted_intervals.append(str(num))
    
    # Final output format: '{val1,val2,...}'
    result_string = "{" + ",".join(formatted_intervals) + "}"
    
    print("The 11 intervals are:")
    print(result_string)

solve_music_intervals()
<<<
{1,1,0.5,1,0.5,1,1,1,1,0.5,1}
>>>