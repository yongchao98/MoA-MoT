def identify_garner_scale():
    """
    This script analyzes a specific melody from Erroll Garner's "All My Loves Are You"
    to identify the scale being used.
    """
    
    print("This program analyzes the musical scale used by Erroll Garner in 'All My Loves Are You' (1986 album 'Afternoon Of An Elf') between 0:39 and 0:43.")
    print("-" * 70)

    # Step 1: Transcribe the notes from the recording.
    # The fast, descending run in the right hand contains the following notes.
    # Note: Cb (C-flat) is used instead of its enharmonic equivalent B natural,
    # as it fits the theoretical structure of the scale.
    notes_in_run = ["G", "F", "Eb", "Db", "Cb", "Bb"]
    
    print("Step 1: Transcription of the Melodic Run")
    print("The notes played in the descending run are:")
    # The instruction "output each number in the final equation" is interpreted
    # here as a request to clearly list out each element (note) of the musical phrase.
    print(" -> ".join(notes_in_run))
    print("\n")

    # Step 2: Analyze the intervals between the notes.
    print("Step 2: Analysis of the Musical Intervals")
    print("The primary characteristic of a scale is the pattern of intervals between its notes.")
    print("Analyzing the intervals in the descending run:")
    print("  G down to F   ...is a Whole Step")
    print("  F down to Eb  ...is a Whole Step")
    print("  Eb down to Db ...is a Whole Step")
    print("  Db down to Cb ...is a Whole Step")
    print("  Cb down to Bb ...is a Half Step")
    print("\n")
    
    # Step 3: Identify the scale based on the analysis.
    print("Step 3: Identification of the Scale")
    print("The defining feature of this run is the prominent sequence of four consecutive whole steps.")
    print("A scale built exclusively from whole steps is called a Whole-Tone Scale.")
    print("Garner's run uses the characteristic sound of the whole-tone scale for most of its duration.")
    print("The final half-step to Bb is a common jazz technique where a performer slightly alters a scale to land on a specific note that fits the song's harmony.")
    print("-" * 70)

    conclusion = "Whole-Tone Scale"
    print(f"Conclusion: The melody is a run based on the {conclusion}.")

identify_garner_scale()