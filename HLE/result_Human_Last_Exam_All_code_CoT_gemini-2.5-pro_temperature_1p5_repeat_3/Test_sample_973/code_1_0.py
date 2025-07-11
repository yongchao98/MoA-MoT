def solve_music_puzzle():
    """
    Solves the music theory puzzle by determining the song's key
    from the chord progression and finding the final melodic note.
    """
    # Using a note list where 'Bb' is the preferred name for A#
    notes = ['C', 'Db', 'D', 'Eb', 'E', 'F', 'Gb', 'G', 'Ab', 'A', 'Bb', 'B']
    
    # A mapping of all note names (including sharps) to their index (0-11)
    note_map = {
        'C': 0, 'B#': 0, 'C#': 1, 'Db': 1, 'D': 2, 'D#': 3, 'Eb': 3,
        'E': 4, 'Fb': 4, 'F': 5, 'E#': 5, 'F#': 6, 'Gb': 6, 'G': 7,
        'G#': 8, 'Ab': 8, 'A': 9, 'A#': 10, 'Bb': 10, 'B': 11, 'Cb': 11
    }

    # The final chord pair in the given progression for "Happy Birthday".
    final_ii_chord = "Cm7"
    final_V_chord = "F7(9)"

    # The root of the ii chord is the first character.
    ii_root_note = final_ii_chord[0]
    
    # A ii chord's root is the second degree of the scale. The tonic (I)
    # is two semitones (a major second) below the ii chord's root.
    
    # Find the numerical index of the ii chord's root.
    ii_root_index = note_map[ii_root_note]

    # To find the tonic's index, we subtract 2. We add 12 and use modulo 12
    # to ensure the result is a valid index between 0 and 11.
    tonic_index = (ii_root_index - 2 + 12) % 12

    # Get the note name of the tonic from our preferred list.
    tonic_note = notes[tonic_index]
    
    print(f"The final phrase of the song is played over the chords '{final_ii_chord}' and '{final_V_chord}'.")
    print("This is a 'ii-V' progression, which points to the tonic (the home note) of the song's key.")
    print(f"The root note of the 'ii' chord ({final_ii_chord}) is '{ii_root_note}'.")
    print(f"The tonic of the key is two semitones below the root of the 'ii' chord.")
    print("We can find the tonic's index with the following equation, using the chromatic scale where C=0:")
    print(f"Equation: (IndexOf({ii_root_note}) - 2 + 12) % 12 = TonicIndex")
    print(f"Calculation: ({ii_root_index} - 2 + 12) % 12 = {tonic_index}")
    print(f"The note at index {tonic_index} is '{tonic_note}'. So, the key of the song is {tonic_note}.")
    print("\nThe standard melody for 'Happy Birthday' always ends on the tonic note.")
    print(f"Therefore, the note sung on the final word, 'you', is {tonic_note}.")

solve_music_puzzle()
<<<Bb>>>