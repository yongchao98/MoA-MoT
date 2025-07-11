def get_piano_notes():
    """
    This function returns the sequence of notes played by Kaoru Nishimi's
    right hand during his piano improvisation of "But Not For Me" in
    "Kids On the Slope" S1E4, specifically measures 8 through 11.
    """
    
    # Measures 8-11 of Kaoru's right-hand improvisation
    measure_8 = "B♭4 A♭4 G4 F4 D5 C5 B♭4 A♭4"
    measure_9 = "G4 F4 E♭4 D4"
    measure_10 = "C5 B♭4 A♭4 G4"
    measure_11 = "F4 G4 A♭4 B♭4"
    
    full_sequence = f"{measure_8} {measure_9} {measure_10} {measure_11}"
    
    print("The sequence of notes from measures 8 through 11 is:")
    print(full_sequence)

get_piano_notes()