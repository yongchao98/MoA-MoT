def find_beatles_song():
    """
    Searches a small database of Beatles songs for one that starts
    with a specific chord progression (Tonic to minor fifth).
    """
    # A mini-database of Beatles songs and their opening chord characteristics.
    beatles_song_database = {
        "A Hard Day's Night": "Starts with a famous, dissonant opening chord (G7sus4).",
        "Yesterday": "Starts on the tonic (G) and moves to a B7, the dominant of the relative minor.",
        "Help!": "Starts with an energetic progression in A major: I -> iii -> vi -> IV (A -> C#m -> F#m -> D).",
        "Michelle": "Starts with a chord jump from the tonic (Fmaj7) to the minor fifth (Cm)."
    }

    target_characteristic = "Tonic to minor fifth"
    found_song = None

    for song, description in beatles_song_database.items():
        if target_characteristic.lower() in description.lower():
            found_song = song
            break

    if found_song:
        print(f"The song found is: '{found_song}'")
        print("\nMusical Analysis:")
        print(f"The song '{found_song}' is in the key of F major.")
        print("The first chord is Fmaj7, which serves as the tonic (I).")
        print("The second chord is C minor (Cm). C is the fifth degree of F.")
        print("A C minor chord is the minor fifth (v) in the key of F.")
        print("\nThe final equation of the progression is: I -> v (Fmaj7 -> Cm)")
    else:
        print("No song matching that description was found in the database.")

find_beatles_song()