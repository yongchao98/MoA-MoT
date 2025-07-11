def calculate_time_signature():
    """
    Calculates the total number of beats in the provided musical measure
    to determine the time signature. We assume a quarter note gets one beat.
    """

    # Note durations in terms of quarter notes (beats)
    eighth_rest = 0.5
    eighth_note = 0.5
    sixteenth_note = 0.25
    quarter_note = 1.0
    half_note = 2.0

    # --- Upper Staff Calculation ---
    # Beat 1: Eighth rest + Eighth note
    upper_beat_1 = eighth_rest + eighth_note
    # Beat 2: Two Sixteenth notes + Eighth note
    upper_beat_2 = sixteenth_note * 2 + eighth_note
    # Beat 3: Two Sixteenth notes + Eighth note
    upper_beat_3 = sixteenth_note * 2 + eighth_note
    # Beat 4: Four Sixteenth notes
    upper_beat_4 = sixteenth_note * 4

    total_upper_beats = upper_beat_1 + upper_beat_2 + upper_beat_3 + upper_beat_4

    print("--- Analysis of the Musical Measure ---")
    print("\nUpper Staff (Right Hand) Beat Calculation:")
    print(f"Beat 1: {eighth_rest} (rest) + {eighth_note} (note) = {upper_beat_1} beat")
    print(f"Beat 2: {sixteenth_note} + {sixteenth_note} + {eighth_note} = {upper_beat_2} beat")
    print(f"Beat 3: {sixteenth_note} + {sixteenth_note} + {eighth_note} = {upper_beat_3} beat")
    print(f"Beat 4: {sixteenth_note} + {sixteenth_note} + {sixteenth_note} + {sixteenth_note} = {upper_beat_4} beat")
    print(f"Total Upper Staff: {upper_beat_1} + {upper_beat_2} + {upper_beat_3} + {upper_beat_4} = {total_upper_beats} beats")

    # --- Lower Staff Calculation ---
    # Beat 1: Quarter note
    lower_beat_1 = quarter_note
    # Beats 2-3: Half note
    lower_beats_2_3 = half_note
    # Beat 4: Quarter note
    lower_beat_4 = quarter_note

    total_lower_beats = lower_beat_1 + lower_beats_2_3 + lower_beat_4

    print("\nLower Staff (Left Hand) Beat Calculation:")
    print(f"Beat 1: {lower_beat_1} (note) = {lower_beat_1} beat")
    print(f"Beats 2-3: {lower_beats_2_3} (note) = {lower_beats_2_3} beats")
    print(f"Beat 4: {lower_beat_4} (note) = {lower_beat_4} beat")
    print(f"Total Lower Staff: {lower_beat_1} + {lower_beats_2_3} + {lower_beat_4} = {total_lower_beats} beats")

    print("\n--- Conclusion ---")
    if total_upper_beats == 4 and total_lower_beats == 4:
        print("Both staves contain 4 quarter-note beats per measure.")
        print("The time signature is 4/4.")
    else:
        print("Could not determine the time signature from the beats.")

calculate_time_signature()