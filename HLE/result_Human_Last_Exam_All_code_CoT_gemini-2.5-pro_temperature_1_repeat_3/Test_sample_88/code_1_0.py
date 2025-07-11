def get_piano_improvisation_notes():
    """
    This function returns the sequence of notes played by Kaoru Nishimi's
    right hand during the improvisation of "But Not For Me" in
    "Kids On the Slope" Season 1 Episode 4, measures 8-11.
    """
    # Notes for Measure 8
    measure_8 = ["G4", "F4", "E♭4", "D4", "C4", "B♭3"]
    
    # Notes for Measure 9
    measure_9 = ["C4", "B♭3", "A♭3", "G3"]
    
    # Notes for Measure 10
    measure_10 = ["G4", "E♭5", "D5", "B♭4", "G4", "E♭4", "D4"]
    
    # Notes for Measure 11
    measure_11 = ["C5", "B♭4", "G4", "E♭4"]
    
    # Combine all measures into a single sequence
    full_sequence = measure_8 + measure_9 + measure_10 + measure_11
    
    # Print the full sequence of notes
    print(" ".join(full_sequence))

get_piano_improvisation_notes()