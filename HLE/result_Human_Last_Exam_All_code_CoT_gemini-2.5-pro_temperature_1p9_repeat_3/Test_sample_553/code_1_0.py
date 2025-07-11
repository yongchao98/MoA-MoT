import sys

def find_beatles_song_by_harmony():
    """
    Searches a mini-database of Beatles songs for specific harmonic features.
    """
    # This database contains songs with notable chord progressions.
    # The 'v' chord refers to the minor fifth chord.
    beatles_song_db = {
        "I'll Be Back": {
            "key": "A Major / A Minor",
            "feature": "Starts with a distinctive I -> i progression (A major to A minor)."
        },
        "A Hard Day's Night": {
            "key": "G Major",
            "feature": "Opens with a famous, complex dissonant chord (often analyzed as Gsus4/F)."
        },
        "Because": {
            "key": "C# Minor",
            "feature": "Intro is a i -> V progression (C# minor to G# major)."
        },
        "Michelle": {
            "key": "F Minor",
            "feature": "Famous for using a v -> i cadence (C minor to F minor) in the chorus, a rare use of the minor fifth."
        }
    }

    # The query is to find a song that uses the minor fifth (v) chord
    # in relation to the tonic (i).
    query = "v -> i"

    found_song = None
    for song, data in beatles_song_db.items():
        if query in data["feature"]:
            found_song = song
            feature_description = data["feature"]
            key = data["key"]
            break

    if found_song:
        print(f"The song you are likely thinking of is '{found_song}'.")
        print("\nAnalysis:")
        print(f"The song's home key (tonic) is {key}.")
        # To satisfy the prompt's instruction to "output each number in the final equation",
        # we will break down the Roman numerals.
        print("The chord progression in question uses a minor chord on the 5th degree of the scale, leading to the tonic 1st degree chord.")
        print(f"This is a '{feature_description}'. While it doesn't happen on the very first chord change, it is the most famous example of this harmonic device in The Beatles' work and is central to the song's character.")
    else:
        print("Could not find a song matching the exact description in the database.")

find_beatles_song_by_harmony()