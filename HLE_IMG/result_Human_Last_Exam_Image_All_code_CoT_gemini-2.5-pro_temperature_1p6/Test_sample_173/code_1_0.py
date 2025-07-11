def find_time_signature():
    """
    Calculates the time signature by summing the note durations in the measure.
    We use the left hand (bass clef) as it has a clear and simple rhythm.
    """
    
    # In music notation, we can assign relative values. Let a quarter note be 1 beat.
    quarter_note_beats = 1.0

    # The left hand contains four notes, each is a quarter note.
    beat_1 = quarter_note_beats
    beat_2 = quarter_note_beats
    beat_3 = quarter_note_beats
    beat_4 = quarter_note_beats

    # Sum the beats to find the total beats in the measure.
    total_beats = beat_1 + beat_2 + beat_3 + beat_4

    print("Analyzing the left hand to find the time signature:")
    print("The measure contains four quarter notes.")
    print("Let's sum the duration of each note, where a quarter note equals 1 beat.")
    print("\nFinal Equation:")
    print(f"{beat_1} + {beat_2} + {beat_3} + {beat_4} = {total_beats} beats")

    # Determine the time signature based on the result.
    numerator = int(total_beats)
    denominator = 4  # The beat is a quarter note.

    print(f"\nThe total is {numerator} beats, and the beat unit is a quarter note.")
    print(f"Therefore, the time signature is {numerator}/{denominator}.")

find_time_signature()