def calculate_time_signature():
    """
    Calculates the time signature of a musical measure by summing the note durations.
    Note values are represented as fractions of a whole note.
    - Whole Note = 1
    - Half Note = 1/2
    - Quarter Note = 1/4
    - Sixteenth Note = 1/16
    """

    # --- Analysis of the Treble Clef (Top Staff) ---
    print("Analyzing the Treble Clef (Top Staff):")
    
    # The top staff consists of four groups of four sixteenth notes.
    num_groups_treble = 4
    notes_per_group_treble = 4
    note_value_treble = 1/16  # Sixteenth note

    total_sixteenth_notes_treble = num_groups_treble * notes_per_group_treble
    total_duration_treble = total_sixteenth_notes_treble * note_value_treble

    # Convert to quarter notes for clarity
    quarter_notes_treble = total_duration_treble / (1/4)
    
    print(f"The top staff has {num_groups_treble} groups of {notes_per_group_treble} sixteenth notes.")
    print(f"Total sixteenth notes = {num_groups_treble} * {notes_per_group_treble} = {total_sixteenth_notes_treble}")
    print(f"This is equivalent to {int(quarter_notes_treble)} quarter notes.\n")

    # --- Analysis of the Bass Clef (Bottom Staff) ---
    print("Analyzing the Bass Clef (Bottom Staff):")
    
    # The bottom staff has two quarter notes and one half note.
    num_quarter_notes_bass = 2
    num_half_notes_bass = 1
    
    quarter_note_value = 1/4
    half_note_value = 1/2

    total_duration_bass = (num_quarter_notes_bass * quarter_note_value) + (num_half_notes_bass * half_note_value)
    
    # Convert to quarter notes
    quarter_notes_bass = total_duration_bass / (1/4)
    
    print(f"The bottom staff has {num_quarter_notes_bass} quarter notes and {num_half_notes_bass} half note.")
    print(f"A half note is equal to 2 quarter notes.")
    print(f"Total quarter notes = {num_quarter_notes_bass} + (1 * 2) = {int(quarter_notes_bass)}")
    
    print("\n--- Conclusion ---")
    print("Both staves contain notes that sum to a total of 4 quarter notes per measure.")
    print("This corresponds to a 4/4 time signature.")
    print("Final Equation: 2 quarter notes + 1 half note = 4 beats")
    print("or")
    print("Final Equation: 16 sixteenth notes = 4 beats")


calculate_time_signature()
<<<C>>>