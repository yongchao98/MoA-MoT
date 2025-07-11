def find_time_signature():
    """
    Calculates the time signature by summing the note durations in the measure.
    We will treat a quarter note as having a value of 1 beat.
    """

    # In 4/4 time, a sixteenth note is worth 0.25 beats (1/4 of a quarter note).
    sixteenth_note_value = 0.25

    # Each of the four groups in the right hand consists of notes and/or rests
    # that add up to the duration of four sixteenth notes.
    group_duration = 4 * sixteenth_note_value  # This is equal to 1 quarter note beat

    # There are four such groups in the measure.
    num_groups = 4

    beat1 = group_duration
    beat2 = group_duration
    beat3 = group_duration
    beat4 = group_duration

    # Calculate the total number of beats in the measure
    total_beats = beat1 + beat2 + beat3 + beat4

    # The time signature's top number is the total beats per measure.
    # The bottom number indicates the note value that gets one beat (4 for a quarter note).
    numerator = int(total_beats)
    denominator = 4

    print("Analyzing the rhythm:")
    print("The measure contains 4 groups. Each group's duration is equivalent to a quarter note (1 beat).")
    print("Equation for total beats:")
    print(f"{int(beat1)} + {int(beat2)} + {int(beat3)} + {int(beat4)} = {int(total_beats)} beats")
    print(f"\nThis means there are {numerator} quarter-note beats per measure.")
    print(f"The time signature is {numerator}/{denominator}.")
    print("This corresponds to answer choice C.")

find_time_signature()