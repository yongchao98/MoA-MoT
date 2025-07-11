def get_kaoru_improvisation():
    """
    This script details the notes from Kaoru Nishimi's piano improvisation
    on "But Not For Me" in "Kids on the Slope," S1E4.
    """

    # The notes are transcribed from measures 8 through 11.
    measure_8 = ["C5", "E♭5", "G5", "A♭5", "G5"]
    measure_9 = ["F5", "G5", "A♭5", "G5", "F5", "E♭5", "D5"]
    measure_10 = ["C5", "E♭5", "G5", "B♭4", "C5", "E♭5", "G5", "B♭4"]
    measure_11 = ["C5", "A4", "F5", "E♭5", "C5", "A4", "G4", "F4"]

    # Combine all measures into a single sequence
    full_sequence = measure_8 + measure_9 + measure_10 + measure_11

    # Print the notes for each individual measure as a step-by-step breakdown
    print("Step 1: Notes for Measure 8")
    for note in measure_8:
        print(note, end=' ')
    print("\n")

    print("Step 2: Notes for Measure 9")
    for note in measure_9:
        print(note, end=' ')
    print("\n")

    print("Step 3: Notes for Measure 10")
    for note in measure_10:
        print(note, end=' ')
    print("\n")
    
    print("Step 4: Notes for Measure 11")
    for note in measure_11:
        print(note, end=' ')
    print("\n")

    # Print the final combined sequence
    print("Final Sequence (Measures 8-11):")
    final_output = " ".join(full_sequence)
    print(final_output)

# Execute the function to get the notes
get_kaoru_improvisation()
<<<C5 E♭5 G5 A♭5 G5 F5 G5 A♭5 G5 F5 E♭5 D5 C5 E♭5 G5 B♭4 C5 E♭5 G5 B♭4 C5 A4 F5 E♭5 C5 A4 G4 F4>>>