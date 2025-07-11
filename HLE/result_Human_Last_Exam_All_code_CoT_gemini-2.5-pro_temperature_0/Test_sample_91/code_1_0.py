def demonstrate_whole_tone_scales():
    """
    This function explains and demonstrates the whole-tone scale,
    which is used by Erroll Garner in the specified musical passage.
    """
    print("Analysis of the musical passage (0:39-0:43) in Erroll Garner's 'All My Loves Are You':")
    print("The right-hand melody in this section is a fast, descending run.")
    print("The sound is characteristic of a WHOLE-TONE SCALE.\n")
    print("A whole-tone scale is built entirely from whole-step intervals.")
    print("This creates a unique, 'dreamy' sound. There are only two possible whole-tone scales:\n")

    # Define the two whole-tone scales
    # Scale 1 starts on C
    whole_tone_scale_1 = ["C", "D", "E", "F#", "G#", "A#"]
    # Scale 2 starts on C#
    whole_tone_scale_2 = ["C#", "D#", "F", "G", "A", "B"]

    print("Whole-Tone Scale 1:")
    print(" -> ".join(whole_tone_scale_1))
    print("\nWhole-Tone Scale 2:")
    print(" -> ".join(whole_tone_scale_2))

    print("\nErroll Garner uses a descending run based on one of these whole-tone scales for his virtuosic flourish.")

# Execute the function to display the analysis
demonstrate_whole_tone_scales()