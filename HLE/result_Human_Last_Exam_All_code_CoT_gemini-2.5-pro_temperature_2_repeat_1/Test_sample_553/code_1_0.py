def find_song_with_specific_chord_jump():
    """
    This script finds a Beatles song from a small, predefined list
    that starts with a specific harmonic movement: a jump from the tonic
    chord (I) to the minor fifth chord (v).
    """
    # A small database of Beatles songs, their key, and their first two chords.
    song_database = {
        "Help!": {"key": "A", "chords": ["A", "C#m"]},
        "A Hard Day's Night": {"key": "G", "chords": ["G11", "C"]},
        "I'll Be Back": {"key": "A", "chords": ["A", "Em"]},
        "Michelle": {"key": "F", "chords": ["C", "Gm7"]},
    }

    print("Searching for a Beatles song with a Tonic (I) to minor fifth (v) chord jump...")

    # --- Music Theory Logic ---
    # The major scale degrees (1-7) for mapping notes.
    scale_map = ['A', 'B', 'C', 'D', 'E', 'F', 'G']

    # This function calculates the note at a specific interval from a root note.
    def get_note_from_interval(root_note, interval_steps):
        start_index = scale_map.index(root_note.upper())
        # The calculation wraps around the scale using the modulo operator.
        target_index = (start_index + interval_steps) % len(scale_map)
        return scale_map[target_index]

    # Iterate through the database to find a match.
    found_song = None
    for song, data in song_database.items():
        key_note = data["key"]
        first_chord = data["chords"][0]
        second_chord = data["chords"][1]

        # 1. The Tonic (I) chord's root note is the key's note itself.
        expected_tonic_note = key_note
        # 2. The fifth (V) chord's root note is 4 steps up the scale map (e.g., A -> E).
        expected_fifth_note = get_note_from_interval(key_note, 4)

        # We are looking for a major Tonic chord and a minor fifth chord.
        expected_tonic_chord = expected_tonic_note
        expected_minor_fifth_chord = expected_fifth_note + "m"

        # Check if the song's chords match the calculated I -> v progression.
        if first_chord == expected_tonic_chord and second_chord == expected_minor_fifth_chord:
            found_song = song
            print("\n--- Match Found! ---")
            print(f"Song: '{found_song}'")
            print(f"Key: {key_note} Major")
            print("\nChord Progression Analysis:")
            print(f"Tonic (I) Chord: {first_chord}")
            print(f"Minor Fifth (v) Chord: {second_chord}")
            print(f"Final Equation: {first_chord} -> {second_chord}")
            break

    if not found_song:
        print("\nNo song in our database matches this specific chord progression.")

# Run the search function.
find_song_with_specific_chord_jump()