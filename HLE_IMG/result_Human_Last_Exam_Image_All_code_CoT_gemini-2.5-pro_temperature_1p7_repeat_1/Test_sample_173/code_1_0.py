def find_time_signature():
    """
    Calculates the time signature by analyzing the note durations in the provided musical measure.
    """

    # --- Note Value Definitions ---
    # We define values relative to a quarter note, which is our 'beat'.
    # A sixteenth note is 1/4 of a quarter note.
    # A whole note is 4 quarter notes.
    sixteenth_note_in_beats = 1.0 / 4.0
    whole_note_in_beats = 4.0

    # --- Analysis of Treble Clef (Top Staff) ---
    # The treble clef contains 16 sixteenth notes.
    num_sixteenth_notes = 16
    total_beats_treble = num_sixteenth_notes * sixteenth_note_in_beats

    # --- Analysis of Bass Clef (Bottom Staff) ---
    # The bass clef contains 1 whole note.
    num_whole_notes = 1
    total_beats_bass = num_whole_notes * whole_note_in_beats

    print("Analyzing the duration of notes to find the time signature.")
    print("The beat unit is assumed to be a quarter note.\n")

    print("--- Treble Clef (Top Staff) Calculation ---")
    print(f"Number of sixteenth notes: {num_sixteenth_notes}")
    print(f"Total beats = {num_sixteenth_notes} * (1/4) = {int(total_beats_treble)} beats\n")

    print("--- Bass Clef (Bottom Staff) Calculation ---")
    print(f"Number of whole notes: {num_whole_notes}")
    print(f"Total beats = {num_whole_notes} * 4 = {int(total_beats_bass)} beats\n")

    # The top number of the time signature is the total beats per measure.
    # The bottom number indicates the note type for one beat (4 for a quarter note).
    top_number = int(total_beats_treble)
    bottom_number = 4

    print("--- Conclusion ---")
    print(f"The measure contains a total of {top_number} quarter-note beats.")
    print("This corresponds to a time signature where the top number is 4 and the bottom number is 4.")
    print(f"\nFinal Equation: The time signature is {top_number} / {bottom_number}")

find_time_signature()