def solve_music_puzzle():
    """
    This script determines the final note of "Happy Birthday" based on its
    initial chord progression.
    """
    
    # The first two chords are Cm7 and F7(9). This is a ii-V progression.
    # The 'V' chord is F7(9), and its root note is 'F'.
    v_chord_root = 'F'
    
    print("The plan is to determine the song's key from its opening chords.")
    print("The first progression is Cm7 -> F7(9), a 'ii-V' progression.")
    print(f"This progression strongly resolves to a 'I' (tonic) chord.")
    print(f"The root of the 'V' chord, F7(9), is '{v_chord_root}'.\n")

    # We can find the tonic by moving up a perfect fourth (5 semitones) from the root of the V chord.
    # Let's define the 12 notes of the chromatic scale.
    chromatic_scale = ['A', 'Bb', 'B', 'C', 'Db', 'D', 'Eb', 'E', 'F', 'Gb', 'G', 'Ab']
    
    # Find the index of the root note.
    root_index = chromatic_scale.index(v_chord_root)
    
    # Calculate the index of the tonic note.
    interval = 5  # A perfect fourth is 5 semitones.
    total_notes = len(chromatic_scale)
    tonic_index = (root_index + interval) % total_notes
    
    # Find the tonic note from the scale.
    tonic_note = chromatic_scale[tonic_index]

    print("The tonic is found by moving up a perfect fourth (5 semitones) from the V-chord's root.")
    print("The calculation for the tonic note's position in the scale is:")
    print(f"(Index of '{v_chord_root}' + Interval) % Total Notes")
    print(f"({root_index} + {interval}) % {total_notes} = {tonic_index}")
    print(f"The note at index {tonic_index} is '{tonic_note}'. Therefore, the song's key is {tonic_note} Major.\n")

    # The melody of "Happy Birthday" always ends on the tonic note.
    final_note = tonic_note
    
    print("The traditional melody of 'Happy Birthday to You' concludes on the tonic note of the key.")
    print(f"Therefore, the note sung on the final word 'you' is:")
    print(final_note)

solve_music_puzzle()
<<<Bb>>>