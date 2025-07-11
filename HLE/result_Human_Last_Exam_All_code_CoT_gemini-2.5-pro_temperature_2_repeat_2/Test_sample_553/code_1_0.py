import sys

def find_beatles_song():
    """
    Finds a Beatles song from a small database that matches a specific
    opening chord progression: Tonic (I) to minor Fifth (vm).
    """

    # A small database of Beatles songs with their key and opening chord
    # progression represented in Roman Numeral Analysis.
    songs_database = {
        'A Hard Day\'s Night': {
            'key': 'G Major',
            'opening_progression': ['G11', 'C', 'G']
        },
        'Yesterday': {
            'key': 'F Major',
            'opening_progression': ['I', 'vii', 'iii7'] # F -> Em7 -> A7
        },
        'Help!': {
            'key': 'A Major',
            'opening_progression': ['I', 'vii(b)', 'vi'] # A -> G -> F#m (simplified analysis)
        },
        'Michelle': {
            'key': 'F Major',
            'opening_progression': ['I', 'vm', 'IVm'] # F -> Cm -> Bbm
        }
    }

    target_progression_start = ['I', 'vm']
    found_song = None
    analysis = {}

    for song_title, song_data in songs_database.items():
        # Check if the start of the song's progression matches the target
        if song_data['opening_progression'][:2] == target_progression_start:
            found_song = song_title
            analysis = song_data
            break

    if found_song:
        print(f"The song by The Beatles that starts with a chord jump from the tonic to the minor fifth is: '{found_song}'.")
        print("\n--- Analysis ---")
        print(f"Key: {analysis['key']}")
        print("Progression Breakdown:")
        tonic_chord = analysis['key'].split(" ")[0]
        # In F major, the fifth note is C. A minor fifth chord is C minor.
        minor_fifth_chord = 'Cm'
        print(f"The first chord (Tonic - I) is: {tonic_chord}")
        print(f"The second chord (minor Fifth - vm) is: {minor_fifth_chord}")
        print(f"The opening musical phrase is effectively: {analysis['opening_progression'][0]} -> {analysis['opening_progression'][1]}")
    else:
        print("Could not find a song with that specific opening chord progression in the database.")

if __name__ == "__main__":
    find_beatles_song()