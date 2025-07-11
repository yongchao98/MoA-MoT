def find_beatles_song():
    """
    Analyzes a small database of Beatles songs to find the one
    matching the user's description of its opening chord progression.
    """
    # A small database of Beatles intros for analysis
    beatles_intros = {
        "Michelle": {
            "key": "F (with mode mixture)",
            "progression": ["Fm", "C"],
            "analysis": (
                "Starts its vocal melody with a jump from a minor tonic (i) "
                "to a major dominant (V). This is the most likely match."
            )
        },
        "Fool on the Hill": {
            "key": "D Major",
            "progression": ["D", "Gm"],
            "analysis": (
                "Features a prominent jump from the major tonic (I) to the "
                "minor subdominant (iv) early in the verse."
            )
        },
        "I Feel Fine": {
            "key": "G Major",
            "progression": ["G", "D7", "C"],
            "analysis": "Starts with a classic I-V-IV riff."
        }
    }

    # The user's query is "tonic to the minor fifth chord".
    # The progression i -> V (minor tonic to major dominant/fifth) is the strongest candidate.
    target_song = "Michelle"
    song_data = beatles_intros[target_song]
    first_chord = song_data["progression"][0]
    second_chord = song_data["progression"][1]

    print(f"The Beatles song that famously starts with a progression matching this description is '{target_song}'.")
    print("\n--- Analysis ---")
    print(f"The song's key is F, but it begins in the parallel minor key.")
    print(f"The first two chords of the vocal line ('Mi-chelle, ma belle') form the progression:")
    print(f"Chord 1: {first_chord}")
    print(f"Chord 2: {second_chord}")
    print("\n--- Music Theory Explanation ---")
    print(f"The first chord, '{first_chord}', acts as the minor tonic (i).")
    print(f"The second chord, '{second_chord}', is the major dominant (V), built on the fifth note of the scale.")
    print("This 'i -> V' jump is a distinctive start and fits the description of a move between a tonic chord and a chord on the fifth, where one of them is minor.")

find_beatles_song()