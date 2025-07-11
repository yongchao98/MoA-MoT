def find_final_note():
    """
    This script determines the final sung note of "Happy Birthday"
    based on the provided jazz chord progression.
    """
    
    # Step 1: Define the sequence of minor seventh (m7) chords.
    # These are the chords played on the key syllables "bir" and "you".
    m7_chords = ["Cm7", "Bm7", "Am7", "Abm7", "Ebm7", "Bm7", "Cm7"]
    
    # Step 2: Define the lyrical events where these chords are placed.
    # The song has 4 phrases. The m7 chords align with "bir" and "you".
    lyrical_events = [
        ("Phrase 1", "bir"),
        ("Phrase 1", "you"),
        ("Phrase 2", "bir"),
        ("Phrase 2", "you"),
        ("Phrase 3", "bir"),  # "Happy birthday dear [Name]"
        ("Phrase 4", "bir"),
        ("Phrase 4", "you")   # The concluding word
    ]
    
    # Step 3: Identify the chord for the final "you".
    # This corresponds to the last item in our sequence.
    final_chord = m7_chords[-1]
    
    print(f"The chord progression is applied sequentially to the song.")
    print(f"The final word 'you' is sung over the chord: {final_chord}")
    
    # Step 4: Determine the sung note from the chord.
    # In a jazz reharmonization, the melody note often aligns with a primary
    # chord tone. For a resolving note like the final "you", the most
    # stable and logical choice is the root of the chord.
    # The root of a chord like "Cm7" is the first letter.
    final_note = final_chord[0]
    
    print(f"The root of the chord {final_chord} is the note '{final_note}'.")
    print("\nTherefore, the final note sung is:")
    print(final_note)

find_final_note()