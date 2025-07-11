def find_final_note():
    """
    Determines the final sung note of "Happy Birthday" based on its chord progression.
    The logic combines music theory concepts regarding key signatures and melody.
    """
    
    # 1. Analyze the chords to find the key.
    # The song starts and ends with the progression Cm7 -> F7(9).
    ii_chord = "Cm7"
    v_chord = "F7(9)"
    
    # This is a ii-V progression. In music theory, a ii-V progression
    # strongly resolves to its 'I' chord, which is the tonic of the key.
    # The root of the 'I' chord is a whole step (2 semitones) below the root of the 'ii' chord.
    ii_chord_root = "C"
    
    # A whole step below 'C' is 'Bb'. So, the key of the song is Bb Major.
    song_key_tonic = "Bb"
    
    # 2. Analyze the melody.
    # The traditional melody of "Happy Birthday to You" always ends on the
    # tonic note of the key. This provides a sense of closure.
    final_melody_note = song_key_tonic

    # 3. Present the conclusion.
    print(f"The plan is to find the key of the song and then identify the final note of the melody.")
    print("-" * 20)
    print(f"Step 1: Determine the song's key from the chords.")
    print(f"The chord progression starts with {ii_chord} -> {v_chord}.")
    print(f"This is a 'ii-V' progression, which points to a tonic key a whole step below the root of the '{ii_chord}'.")
    print(f"The root of {ii_chord} is {ii_chord_root}.")
    print(f"Calculation: Key = Root({ii_chord}) - 1 whole step")
    print(f"Result: The key is {song_key_tonic} Major.")
    print("-" * 20)
    print(f"Step 2: Determine the final note of the melody.")
    print(f"The song 'Happy Birthday to You' traditionally ends on the tonic note of its key.")
    print(f"Final Note = Tonic of the Key")
    print("-" * 20)
    print(f"Conclusion: Since the key is {song_key_tonic} Major, the final note sung on the last word 'you' is the tonic.")
    print(f"\nThe final note is: {final_melody_note}")

# Run the analysis to find the final note.
find_final_note()