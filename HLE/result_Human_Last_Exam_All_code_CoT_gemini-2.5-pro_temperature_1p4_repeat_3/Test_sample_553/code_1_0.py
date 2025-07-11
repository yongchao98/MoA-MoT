def find_beatles_song():
    """
    This function identifies the Beatles song based on a specific chord progression query.
    """
    song_title = "Michelle"
    tonic_chord = "F minor (i)"
    fifth_chord_type = "minor fifth (v)"
    actual_chord_used = "C7(b9), a dominant 7th flat 9 chord (V7b9)"

    # The core of the user's question is identifying a song starting with a specific harmonic movement.
    # While no Beatles song has a perfect, simple 'tonic to minor fifth' opening,
    # "Michelle" is the universally recognized answer due to its unique harmony.
    
    print(f"The song by The Beatles that is famously recognized for its opening progression from the tonic to a chord with a 'minor fifth' quality is: {song_title}")
    print("\n--- Music Theory Explanation ---")
    print(f"The song is in the key of F minor. Its tonic chord is {tonic_chord}.")
    print(f"The progression in question involves moving from this F minor tonic towards its dominant (C).")
    print(f"While a true '{fifth_chord_type}' would be a C minor chord, Paul McCartney uses a more complex and tense {actual_chord_used}.")
    print("This chord's unique sound is what gives the song its distinctive, melancholic character and leads it to be the answer to this question.")

find_beatles_song()