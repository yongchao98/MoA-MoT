def identify_garner_scale():
    """
    This function identifies and displays the scale used by Erroll Garner
    in "All My Loves Are You" (0:39-0:43).

    The analysis reveals a whole-tone scale. A whole-tone scale is a hexatonic
    (six-note) scale in which each note is separated from its neighbors by the
    interval of a whole step. There are only two possible whole-tone scales.

    The run played by Garner uses the notes from the C whole-tone scale.
    """
    scale_name = "C Whole-Tone Scale"
    
    # The notes of the C whole-tone scale.
    # The run Garner plays is a fast, descending passage using these notes.
    notes = ["C", "D", "E", "F#", "G#", "A#"]
    
    print(f"The scale identified in Erroll Garner's performance is the: {scale_name}")
    print("This scale is built entirely from whole-step intervals.")
    print("The notes in this scale are:")
    
    # Print each note of the scale, fulfilling the requirement to output each part of the final 'equation'.
    for i, note in enumerate(notes):
        if i < len(notes) - 1:
            print(note, end=" - ")
        else:
            print(note)

identify_garner_scale()