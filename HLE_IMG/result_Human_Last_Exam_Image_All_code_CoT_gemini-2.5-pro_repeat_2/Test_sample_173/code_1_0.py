def calculate_time_signature():
    """
    Calculates the time signature of a musical measure by summing the duration of its notes.
    """
    # Let's define the duration of notes relative to a quarter note beat.
    # A quarter note gets 1 beat.
    sixteenth_note_duration = 0.25  # A sixteenth note is 1/4 of a quarter note.
    sixteenth_rest_duration = 0.25

    # Count the notes and rests in the top staff of the measure.
    num_sixteenth_notes = 15
    num_sixteenth_rests = 1

    # Calculate the total duration in terms of quarter note beats.
    total_beats_from_notes = num_sixteenth_notes * sixteenth_note_duration
    total_beats_from_rests = num_sixteenth_rests * sixteenth_rest_duration
    total_beats = total_beats_from_notes + total_beats_from_rests

    print("Analyzing the top staff to determine the time signature.")
    print(f"The measure contains {num_sixteenth_notes} sixteenth notes and {num_sixteenth_rests} sixteenth rest.")
    print("We will calculate the total number of quarter-note beats.\n")
    print("Equation:")
    print(f"Total Beats = (Number of sixteenth notes * Duration) + (Number of sixteenth rests * Duration)")
    print(f"Total Beats = ({num_sixteenth_notes} * {sixteenth_note_duration}) + ({num_sixteenth_rests} * {sixteenth_rest_duration})")
    print(f"Total Beats = {total_beats_from_notes} + {total_beats_from_rests}")
    print(f"Total Beats = {total_beats}\n")
    
    # The top number of the time signature is the total number of beats.
    # The bottom number represents the note value that gets one beat (4 for a quarter note).
    top_number = int(total_beats)
    bottom_number = 4

    print(f"The calculation shows there are {top_number} quarter-note beats per measure.")
    print(f"Therefore, the time signature is {top_number}/{bottom_number}.")
    print("\nThis is confirmed by the single whole note in the bottom staff, which also equals 4 quarter notes.")

calculate_time_signature()