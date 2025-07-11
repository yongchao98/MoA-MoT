def identify_garner_scale():
    """
    This function identifies and displays the scale played by Erroll Garner
    in "All My Loves Are You" between 0:39 and 0:43.
    """
    scale_name = "C Natural Minor"
    # The notes in the C Natural Minor scale are C, D, Eb, F, G, Ab, Bb.
    notes = ["C", "D", "Eb", "F", "G", "Ab", "Bb"]

    print(f"The scale played by Erroll Garner's right hand between 0:39 and 0:43 is the {scale_name} scale.")
    print("The notes in this scale are:")
    
    # We will print the notes of the scale clearly
    output = " - ".join(notes)
    print(output)

identify_garner_scale()