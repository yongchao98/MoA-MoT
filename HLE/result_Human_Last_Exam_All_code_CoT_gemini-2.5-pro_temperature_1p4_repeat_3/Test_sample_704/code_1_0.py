def find_newtonian_note():
    """
    This function determines the musical note corresponding to a specific color range
    on Isaac Newton's color circle.
    """
    # According to Isaac Newton's Opticks, he mapped the notes of the
    # diatonic scale to the boundaries between the seven colors of the spectrum.
    newton_color_boundaries_to_notes = {
        ('Red', 'Orange'): 'D',
        ('Orange', 'Yellow'): 'E',
        ('Yellow', 'Green'): 'F',
        ('Green', 'Blue'): 'G',
        ('Blue', 'Indigo'): 'A',
        ('Indigo', 'Violet'): 'B',
        ('Violet', 'Red'): 'C'
    }

    # The problem asks for the note in the range between the combination
    # of yellow and blue (which is green) and blue itself.
    start_color_of_range = "Green"
    end_color_of_range = "Blue"

    # We look for the note at the boundary between Green and Blue.
    boundary_key = (start_color_of_range, end_color_of_range)

    # Retrieve the note from the mapping.
    musical_note = newton_color_boundaries_to_notes[boundary_key]

    print("Newton's Color Circle maps colors to musical notes.")
    print("The starting color of the range is the combination of Yellow and Blue, which is Green.")
    print("The ending color of the range is Blue.")
    print(f"The musical note at the boundary between {start_color_of_range} and {end_color_of_range} is:")
    print(musical_note)

find_newtonian_note()