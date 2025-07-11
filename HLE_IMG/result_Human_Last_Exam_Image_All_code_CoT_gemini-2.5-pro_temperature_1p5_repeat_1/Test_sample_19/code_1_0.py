def moonligh_sonata_analysis():
    """
    This function provides a three-part music theory analysis of the provided
    excerpt from Beethoven's Moonlight Sonata, 1st movement.
    """

    # Answer to part 1
    answer1 = ("1. To what key the music modulates in measure 11 and the beginning of measure 12?\n\n"
               "The music modulates back to the home key of C-sharp minor. "
               "Measure 11 contains a highly chromatic passage that creates tension, "
               "which is then resolved with the arrival of a clear C-sharp minor tonic chord "
               "on the first beat of measure 12.")

    # Answer to part 2
    answer2 = ("2. Given the 4 sharps key environment, what is the connection and / or justification for this "
               "particular modulation?\n\n"
               "The key signature of four sharps can be either E major or C-sharp minor. The piece's home key is C-sharp minor. "
               "E major is the relative major of C-sharp minor. In the measures preceding this passage, the music briefly tonicizes "
               "E major. The 'modulation' in measure 11 is therefore a return to the tonic (i). Returning to the home key, especially "
               "from its closely-related relative major, is a fundamental and logical harmonic progression that reinforces the tonal center.")

    # Answer to part 3
    answer3 = ("3. What would be the complete and accurate roman numeral marking for the first beat of measure 11?\n\n"
               "The complete Roman numeral marking is vii°65/iv (vii dim 6-5 of iv).\n"
               "Here is the justification:\n"
               " - The chord contains the notes B-natural, D-natural, and F-natural.\n"
               " - Its function is a secondary (or applied) leading-tone chord that points to the subdominant chord (iv), which is F-sharp minor.\n"
               " - The proper leading-tone seventh chord of F-sharp minor is E-sharp diminished seventh (notes: E#, G#, B, D).\n"
               " - Beethoven writes the chord using an enharmonic spelling (F-natural instead of E-sharp) and omits the G-sharp.\n"
               " - The bass note is B, which is the third of the chord's root (E#), making it a first inversion chord. The figured bass for a first inversion seventh chord is 6/5 (often written as 65).")
    
    # Print all answers
    print(answer1)
    print("\n----------------------------------------\n")
    print(answer2)
    print("\n----------------------------------------\n")
    print(answer3)

moonligh_sonata_analysis()