def calculate_time_signature():
    """
    Calculates the time signature by summing the note durations in the provided musical measure.
    Note values are represented in terms of quarter notes.
    1 sixteenth note = 0.25 quarter notes
    1 half note = 2 quarter notes
    """
    # Analysis of the top stave (right hand)
    print("Analyzing the top stave (right hand):")
    num_sixteenth_notes_per_group = 4
    num_groups = 4
    total_sixteenth_notes = num_sixteenth_notes_per_group * num_groups
    
    # There are 4 sixteenth notes in a quarter note.
    quarters_in_a_sixteenth = 1 / 4
    total_quarter_notes_rh = total_sixteenth_notes * quarters_in_a_sixteenth
    
    print(f"The top stave has {num_groups} groups of {num_sixteenth_notes_per_group} sixteenth notes.")
    print(f"Calculation: {num_groups} * {num_sixteenth_notes_per_group} * (1/4) = {int(total_quarter_notes_rh)} quarter notes.")

    # Analysis of the bottom stave (left hand)
    print("\nAnalyzing the bottom stave (left hand):")
    num_half_notes = 2
    quarters_in_a_half_note = 2
    total_quarter_notes_lh = num_half_notes * quarters_in_a_half_note
    
    print(f"The bottom stave has {num_half_notes} half notes.")
    print(f"Calculation: {num_half_notes} * {quarters_in_a_half_note} = {total_quarter_notes_lh} quarter notes.")

    print(f"\nConclusion: The measure contains {int(total_quarter_notes_rh)} quarter notes, so the time signature is 4/4.")

calculate_time_signature()