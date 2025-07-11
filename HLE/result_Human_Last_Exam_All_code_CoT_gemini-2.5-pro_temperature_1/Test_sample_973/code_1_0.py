def solve_music_puzzle():
    """
    Determines the final sung note of "Happy Birthday" based on harmonic context.
    """
    # Step 1: Define musical constants and given information.
    
    # The C Major Scale, mapping scale degrees to note names.
    # The tonic (or root) is at degree 1.
    c_major_scale = {1: "C", 2: "D", 3: "E", 4: "F", 5: "G", 6: "A", 7: "B"}
    
    # In music theory, the "Happy Birthday" melody starts on the 5th degree
    # (the dominant) and ends on the 1st degree (the tonic).
    starting_melody_degree = 5
    final_melody_degree = 1

    # The first chord given in the jazz arrangement.
    first_chord_name = "Cm7"
    first_chord_notes = {"C", "Eb", "G", "Bb"}

    print("Step 1: Determine the musical key of the song.")
    
    # Find the note name corresponding to the starting melody degree.
    starting_note_name = c_major_scale[starting_melody_degree]
    
    print(f"The standard melody for 'Happy Birthday' starts on the 5th degree of the scale.")
    print(f"In the key of C, the 5th degree is the note '{starting_note_name}'.")
    print(f"The first chord played in the arrangement is {first_chord_name}, which contains the notes {sorted(list(first_chord_notes))}.")

    # Check if the starting melody note is part of the first chord.
    if starting_note_name in first_chord_notes:
        print(f"The starting melody note, '{starting_note_name}', is found within the first chord, {first_chord_name}.")
        print("This strong harmonic relationship confirms that the song is in the key of C.")
    else:
        # This case is not expected here, but included for completeness.
        print("Analysis inconclusive: The starting melody note does not fit the first chord.")
        return

    print("\nStep 2: Identify the note for the final word 'you'.")
    print("The 'Happy Birthday' melody traditionally ends on the tonic (1st degree) of the key.")
    
    # The final note is the tonic of the established key.
    final_note = c_major_scale[final_melody_degree]

    print("\nFinal Equation:")
    print(f"Final Sung Note = Note at Scale Degree {final_melody_degree} in the Key of C")
    print(f"Final Sung Note = {final_note}")

solve_music_puzzle()