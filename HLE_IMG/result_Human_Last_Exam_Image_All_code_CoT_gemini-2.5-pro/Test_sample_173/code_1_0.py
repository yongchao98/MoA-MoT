def calculate_time_signature():
    """
    Calculates the time signature by summing the note durations in the bass clef.
    In this measure, the bass clef contains four quarter notes.
    We'll represent a quarter note as 1 beat.
    """
    # Number of quarter notes in the bass clef
    num_quarter_notes = 4
    
    # Value of a single quarter note in beats
    quarter_note_value = 1
    
    # List of the notes' beat values
    notes_in_measure = [quarter_note_value] * num_quarter_notes
    
    # Calculate the total number of beats
    total_beats = sum(notes_in_measure)
    
    # Print the explanation and the calculation
    print("To find the time signature, we sum the duration of the notes in the measure.")
    print("The bass clef contains 4 quarter notes.")
    print("If a quarter note gets one beat, the calculation is:")
    
    # Create the equation string
    equation = " + ".join(map(str, notes_in_measure))
    
    print(f"{equation} = {total_beats} beats")
    print(f"This corresponds to a 4/4 time signature.")

calculate_time_signature()