def find_final_note():
    """
    Determines the final sung note of "Happy Birthday" based on the
    provided chord progression.
    """
    # 1. The problem provides the chord progression. The final pair given for
    # the last "happy birthday" is Cm7 followed by F7(9).
    final_ii_chord = "Cm7"
    final_V_chord = "F7(9)"

    # 2. This progression, a ii-V, establishes the musical key. The V chord
    # (F7(9)) has a root note of 'F'. The key's tonic is a perfect fifth
    # below the root of the V chord.
    # On a musical keyboard, a perfect fifth is a distance of 7 semitones.
    # The note 7 semitones below 'F' is 'Bb'.
    tonic_note = "Bb"

    # 3. The standard melody of "Happy Birthday" concludes on the tonic note
    # of the key. The final word "you" is sung on this note.
    final_sung_note = tonic_note

    # Print the step-by-step reasoning
    print("Step 1: The final harmonic progression for 'happy birthday' is identified.")
    print(f"The chords are {final_ii_chord} -> {final_V_chord}.")
    print("\nStep 2: The musical key is determined from this progression.")
    print(f"The chord {final_V_chord} (a V chord with root F, a dominant 7th, and a 9th) resolves to the tonic key of {tonic_note}.")
    
    print("\nStep 3: The final note of the melody is identified.")
    print("The song 'Happy Birthday to You' traditionally ends on the tonic note for the final word 'you'.")

    print("\nConclusion:")
    print("Since the key is {} and the song ends on the tonic, the final note sung is {}.".format(tonic_note, final_sung_note))

# Execute the logic and print the result.
find_final_note()