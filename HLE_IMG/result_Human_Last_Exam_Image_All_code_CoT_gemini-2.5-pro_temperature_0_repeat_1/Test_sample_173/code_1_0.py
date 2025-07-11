def calculate_time_signature():
    """
    Calculates the time signature by analyzing the note durations in the provided musical measure.
    """
    # Note durations relative to a whole note.
    note_values = {
        "whole": 1.0,
        "quarter": 0.25,
        "sixteenth": 0.0625
    }

    # --- Right Hand Analysis (Top Staff) ---
    num_sixteenth_notes_rh = 16
    # Total duration of the right hand relative to a whole note.
    total_duration_rh = num_sixteenth_notes_rh * note_values["sixteenth"]
    # Convert the total duration to the number of quarter notes.
    num_quarter_notes_rh = total_duration_rh / note_values["quarter"]

    print("--- Analysis of the Measure ---")
    print("\nRight Hand (Top Staff):")
    print(f"The part contains {num_sixteenth_notes_rh} sixteenth notes.")
    print(f"Calculation: {num_sixteenth_notes_rh} sixteenth notes * (1/16) = {int(num_quarter_notes_rh)} quarter notes.")

    # --- Left Hand Analysis (Bottom Staff) ---
    num_whole_notes_lh = 1
    # Total duration of the left hand relative to a whole note.
    total_duration_lh = num_whole_notes_lh * note_values["whole"]
    # Convert the total duration to the number of quarter notes.
    num_quarter_notes_lh = total_duration_lh / note_values["quarter"]

    print("\nLeft Hand (Bottom Staff):")
    print(f"The part contains {num_whole_notes_lh} whole note.")
    print(f"Calculation: {num_whole_notes_lh} whole note = {int(num_quarter_notes_lh)} quarter notes.")

    # --- Conclusion ---
    print("\n--- Conclusion ---")
    if num_quarter_notes_rh == 4 and num_quarter_notes_lh == 4:
        print("Each measure contains the equivalent of 4 quarter notes.")
        print("This corresponds to a 4/4 time signature.")
    else:
        print("Could not determine the time signature from the notes.")

calculate_time_signature()