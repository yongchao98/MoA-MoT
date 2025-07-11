def find_happy_birthday_note():
    """
    This function solves for the final note of "Happy Birthday" based on the provided
    chord progression. It uses music theory principles to determine the song's key.
    """

    # The problem states the final chord pair is Cm7 -> F7(9).
    # This is a ii-V progression. The 'V' chord is F7. We only need its root note.
    final_V_chord_root = "F"

    # Define the 12 notes of the chromatic scale to perform calculations.
    # We use flats because the key of Bb contains flats.
    chromatic_scale = ["C", "Db", "D", "Eb", "E", "F", "Gb", "G", "Ab", "A", "Bb", "B"]

    # Find the index of the V-chord's root note.
    v_chord_index = chromatic_scale.index(final_V_chord_root)

    # The tonic of a key is a perfect fourth (5 semitones) above the root of its V chord.
    # We can calculate its index using modular arithmetic.
    semitones_in_perfect_fourth = 5
    number_of_notes = 12
    tonic_index = (v_chord_index + semitones_in_perfect_fourth) % number_of_notes

    # Find the note name for the tonic.
    tonic_note = chromatic_scale[tonic_index]

    print("To find the final note, we must first establish the song's key.")
    print("The final chord progression given is Cm7 -> F7(9), which is a 'ii-V' progression.")
    print(f"The V-chord's root is '{final_V_chord_root}'. This chord strongly resolves to the tonic of the key.")
    print(f"The tonic note is 5 semitones (a perfect fourth) higher than the root of the V-chord.")
    print("\nLet's calculate the tonic note:")
    print(f"The index of '{final_V_chord_root}' in the chromatic scale is {v_chord_index}.")
    print(f"The equation to find the tonic's index is: ({v_chord_index} + {semitones_in_perfect_fourth}) % {number_of_notes} = {tonic_index}")
    
    print(f"\nThe note at index {tonic_index} is '{tonic_note}'. So, the key is {tonic_note} Major.")
    print("The melody of 'Happy Birthday' always ends on the tonic note.")
    print("\nTherefore, the note sung on the concluding word 'you' is:")
    print(tonic_note)

find_happy_birthday_note()
<<<Bb>>>