def find_final_note():
    """
    This function determines the final sung note of "Happy Birthday" based on
    the provided chord progression and rules.
    """
    # The chord progression is structured in pairs.
    # Each pair corresponds to a lyrical phrase.
    progression_pairs = [
        ("Cm7", "F7(9)"),   # Phrase 1: "Happy birthday"
        ("Bm7", "E7(9)"),   # Phrase 2: "to you"
        ("Am7", "D7(9)"),   # Phrase 3: "Happy birthday"
        ("Abm7", "Db7(9)"), # Phrase 4: "dear [Name]"
        ("Ebm7", "Ab7(9)"), # Phrase 5: "Happy birthday"
        ("Bm7", "E7(9)")    # Phrase 6: "to you" (FINAL)
    ]

    # The final "you" is sung during the 6th phrase of the song.
    # In programming, this corresponds to index 5.
    final_phrase_index = 5
    print(f"Step 1: The final phrase corresponds to the chord pair at index {final_phrase_index} in the progression.")

    # Get the chord pair for the final phrase.
    final_chord_pair = progression_pairs[final_phrase_index]
    print(f"Step 2: The identified chord pair is ({final_chord_pair[0]}, {final_chord_pair[1]}).")

    # The prompt states the minor 7th chord (the first in the pair) is sung on "you".
    # This corresponds to index 0 of the pair.
    chord_on_you_index = 0
    final_chord = final_chord_pair[chord_on_you_index]
    print(f"Step 3: The chord for the word 'you' is the one at index {chord_on_you_index} of the pair, which is {final_chord}.")

    # The final sung note is typically the root of the chord.
    # The root is the first character of the chord name, at index 0.
    root_note_index = 0
    final_note = final_chord[root_note_index]
    print(f"Step 4: The root note of {final_chord} is the character at index {root_note_index}, which is '{final_note}'.")

    # The final "equation" combines these steps.
    print("\n--- Final Calculation ---")
    print(f"Final Note = progression_pairs[{final_phrase_index}][{chord_on_you_index}][{root_note_index}]")
    print(f"The numbers used in this final equation are {final_phrase_index}, {chord_on_you_index}, and {root_note_index}.")
    print(f"The resulting final note is: {final_note}")

find_final_note()
<<<B>>>