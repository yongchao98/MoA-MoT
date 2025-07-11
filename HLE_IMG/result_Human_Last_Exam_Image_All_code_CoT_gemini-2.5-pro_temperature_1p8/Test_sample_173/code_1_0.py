def find_time_signature():
    """
    Calculates the time signature by summing the note durations in the measure.
    We will analyze the simpler bottom staff (left hand).
    """

    # In music notation, note durations are fractions of a whole note.
    # A quarter note has a value of 1/4.
    quarter_note_value = 1/4
    note_type = "1/4"

    # The bottom staff clearly contains three notes. By their shape and context,
    # each is a quarter note.
    number_of_notes = 3

    # Calculate the total duration of the measure.
    total_duration = number_of_notes * quarter_note_value

    print("To find the time signature, we sum the duration of all notes in the measure.")
    print(f"The bottom staff contains {number_of_notes} notes, each identified as a quarter note ({note_type}).")
    print("\nThe calculation is as follows:")

    # Create the equation string
    equation_parts = [note_type] * number_of_notes
    equation_str = " + ".join(equation_parts)

    print(f"{equation_str} = {total_duration}")

    print("\nA total duration of 3/4 per measure corresponds to the time signature 3/4.")
    print("This means there are 3 beats in the measure, and a quarter note gets one beat.")


find_time_signature()