def find_final_note():
    """
    This function solves the problem by applying music theory principles.
    1. It identifies the musical key from the ii-V chord progression.
    2. It determines the tonic note based on this key.
    3. It concludes the final note based on the standard melody of "Happy Birthday".
    """

    # Using a chromatic scale with flats, where each note has a numeric index.
    # We start with A=0 for this example.
    # A, Bb, B,  C, Db, D, Eb, E, F, Gb, G, Ab
    # 0, 1,  2,  3, 4,  5, 6,  7, 8, 9,  10, 11
    chromatic_scale = ['A', 'Bb', 'B', 'C', 'Db', 'D', 'Eb', 'E', 'F', 'Gb', 'G', 'Ab']

    # The chord progression Cm7 -> F7 is a ii-V progression.
    # The root of the ii chord (Cm7) is C.
    ii_chord_root = 'C'

    # Find the numerical index of the root of the ii chord.
    ii_index = chromatic_scale.index(ii_chord_root)

    # The tonic (I) of the key is 2 semitones (a whole step) below the ii chord's root.
    # We can calculate its index. We use the modulo operator (%) to wrap around the scale if needed.
    tonic_index = (ii_index - 2) % 12

    # Let's display the equation used to find the tonic's position in the scale.
    # We output each number in the equation as requested.
    num1 = ii_index
    num2 = 2
    result = tonic_index
    print(f"The root of the 'ii' chord, '{ii_chord_root}', is at index {num1} in our scale.")
    print(f"The tonic note's index is 2 positions lower.")
    print(f"Final Equation: {num1} - {num2} = {result} (modulo 12)")


    # Find the note name for the tonic from its calculated index.
    tonic_note = chromatic_scale[tonic_index]

    print(f"\nThis calculation shows the key of the song is {tonic_note} Major.")

    # The melody of "Happy Birthday" traditionally ends on the tonic note of the key.
    print(f"The song 'Happy Birthday' melodically resolves to the tonic on the final word, 'you'.")
    print(f"Therefore, the concluding note sung is: {tonic_note}")

find_final_note()
<<<Bb>>>