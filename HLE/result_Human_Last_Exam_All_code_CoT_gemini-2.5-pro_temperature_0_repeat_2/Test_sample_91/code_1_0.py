def identify_garner_scale():
    """
    This function identifies and describes the scale used by Erroll Garner
    in "All My Loves Are You" (0:39-0:43).
    """
    # Analysis of the musical passage reveals a descending run over a G dominant chord.
    # The notes used are characteristic of a specific octatonic scale.
    scale_name = "Half-Whole Diminished Scale"
    root_note = "G"

    # The Half-Whole Diminished scale is an 8-note scale built by
    # alternating half steps and whole steps.
    # Formula starting from the root: Half, Whole, Half, Whole, Half, Whole, Half, Whole
    # G -> Ab (Half)
    # Ab -> Bb (Whole)
    # Bb -> B (Half)
    # B -> C# (Whole)
    # C# -> D (Half)
    # D -> E (Whole)
    # E -> F (Half)
    # F -> G (Whole)
    notes_in_scale = ["G", "Ab", "Bb", "B", "C#", "D", "E", "F"]

    # The run Garner plays is a descending pattern using notes from this scale.
    # The transcribed run is approximately: G, F, E, D, C#, B, Bb (or A#).
    garner_run = ["G", "F", "E", "D", "C#", "B", "Bb"]

    print(f"In the song 'All My Loves Are You' (0:39-0:43), Erroll Garner plays a run based on the {root_note} {scale_name}.")
    print("\nThis scale is frequently used in jazz over dominant 7th chords.")
    print("The notes in the full scale are:")
    # As requested, printing each note of the scale "equation".
    print(' '.join(notes_in_scale))

    print("\nThe specific descending phrase played by Garner consists of the following notes from that scale:")
    print(' '.join(garner_run))

if __name__ == '__main__':
    identify_garner_scale()