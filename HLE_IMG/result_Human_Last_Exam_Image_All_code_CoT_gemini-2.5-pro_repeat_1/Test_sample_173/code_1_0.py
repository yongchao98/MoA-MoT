def find_time_signature():
    """
    Calculates the time signature by analyzing the note durations in the provided musical measure.
    """
    print("To find the time signature, we sum the rhythmic values of the notes in the measure.")
    print("Let's define a quarter note as 1 beat.\n")

    # --- Analysis of the Bass Clef (Bottom Staff) ---
    print("Analyzing the Bass Clef (bottom staff):")
    num_quarter_notes = 4
    beats_per_quarter_note = 1
    total_beats_bass = num_quarter_notes * beats_per_quarter_note
    print(f"The staff contains {num_quarter_notes} quarter notes.")
    print(f"Equation: {num_quarter_notes} (notes) * {beats_per_quarter_note} (beat/note) = {total_beats_bass} beats.")
    print("-" * 30)

    # --- Analysis of the Treble Clef (Top Staff) ---
    print("Analyzing the Treble Clef (top staff):")
    num_sixteenth_notes = 16
    beats_per_sixteenth_note = 1/4  # A sixteenth note is 1/4 of a beat
    total_beats_treble = num_sixteenth_notes * beats_per_sixteenth_note
    print(f"The staff contains {num_sixteenth_notes} sixteenth notes.")
    # For clearer output, we show the division
    print(f"Since 4 sixteenth notes equal 1 beat, we calculate:")
    print(f"Equation: {num_sixteenth_notes} (notes) / 4 (notes/beat) = {int(total_beats_treble)} beats.")
    print("-" * 30)

    # --- Final Conclusion ---
    print(f"Both staves add up to {int(total_beats_bass)} beats.")
    print("This means the measure has 4 beats, and the quarter note gets the beat.")
    print("Therefore, the time signature is 4/4.")

find_time_signature()