def calculate_time_signature():
    """
    Calculates the time signature by summing the note durations in the measure.
    Note values are relative to a quarter note (1.0).
    """
    # Note value definitions
    quarter_note = 1.0
    sixteenth_note = 0.25

    print("--- Analysis of the Musical Measure ---")

    # --- Right Hand (Treble Clef) Analysis ---
    print("\n1. Analyzing the Right Hand (Treble Clef):")
    # The right hand is composed of three groups of four sixteenth notes.
    rh_group_1_value = 4 * sixteenth_note
    rh_group_2_value = 4 * sixteenth_note
    rh_group_3_value = 4 * sixteenth_note
    rh_total_duration = rh_group_1_value + rh_group_2_value + rh_group_3_value

    print("The right hand contains three groups of four sixteenth notes.")
    print(f"The total duration is calculated by summing the value of each group (where 4 sixteenth notes = {4 * sixteenth_note} quarter note).")
    print(f"Equation: {rh_group_1_value} + {rh_group_2_value} + {rh_group_3_value} = {rh_total_duration} quarter notes.")

    # --- Left Hand (Bass Clef) Analysis ---
    print("\n2. Analyzing the Left Hand (Bass Clef):")
    # The left hand is composed of three quarter notes.
    lh_note_1_value = quarter_note
    lh_note_2_value = quarter_note
    lh_note_3_value = quarter_note
    lh_total_duration = lh_note_1_value + lh_note_2_value + lh_note_3_value

    print("The left hand contains three quarter notes.")
    print("The total duration is calculated by summing the value of each note.")
    print(f"Equation: {lh_note_1_value} + {lh_note_2_value} + {lh_note_3_value} = {lh_total_duration} quarter notes.")

    # --- Conclusion ---
    print("\n3. Conclusion:")
    if rh_total_duration == lh_total_duration:
        beats_per_measure = int(rh_total_duration)
        beat_unit = 4  # Since the count is in quarter notes
        print(f"Both hands have a total duration equivalent to {beats_per_measure} quarter notes.")
        print(f"This indicates a time signature of {beats_per_measure}/{beat_unit}.")
    else:
        print("Error: The durations in the hands do not match.")

calculate_time_signature()