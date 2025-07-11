def find_final_note():
    """
    This function determines the final note of "Happy Birthday" based on the provided chord progression.
    
    1. The final sung word is "you".
    2. The chord played over the final "you" is Bm7.
    3. The notes in a Bm7 chord are B, D, F#, and A.
    4. The melody of "Happy Birthday" resolves to the most stable note, the tonic.
    5. In this arrangement, the most logical and musically consonant final note that aligns with the Bm7 chord is the root of the chord itself.
    """
    
    final_chord = "Bm7"
    final_chord_root = "B"
    final_word = "you"
    
    print(f"The final word sung is '{final_word}'.")
    print(f"The chord played over this word is {final_chord}.")
    print(f"For the melody to resolve naturally, it should land on the most stable note of the chord.")
    print(f"The most stable note in the {final_chord} chord is its root.")
    print(f"Therefore, the final note sung is: {final_chord_root}")

find_final_note()
<<<B>>>