def find_final_note():
    """
    This function determines the final sung note of "Happy Birthday to You"
    by analyzing the standard melody of its final phrase.
    """
    
    # The question asks for the sung note, which is determined by the song's melody,
    # not the underlying chords. The standard melody of "Happy Birthday" is consistent.
    # The notes for the final phrase, "Happy birthday to you," are F, F, E, C, D, C.

    # We can represent this melody as a sequence of notes.
    final_phrase_melody = ['F', 'F', 'E', 'C', 'D', 'C']
    
    # We can also represent this as a 'melodic equation' using numbers where C=1, D=2, E=3, F=4.
    melodic_equation = [4, 4, 3, 1, 2, 1]

    # The final word, "you", corresponds to the last note in this sequence.
    final_note = final_phrase_melody[-1]
    
    print("The melodic sequence for the final phrase 'Happy birthday to you' can be represented as the 'equation':")
    # As requested, output each number in the final equation
    print(f"{melodic_equation[0]}, {melodic_equation[1]}, {melodic_equation[2]}, {melodic_equation[3]}, {melodic_equation[4]}, {melodic_equation[5]}")

    print("\nThe concluding word, 'you', is sung on the final note of this sequence.")
    print(f"The final note is: {final_note}")

find_final_note()