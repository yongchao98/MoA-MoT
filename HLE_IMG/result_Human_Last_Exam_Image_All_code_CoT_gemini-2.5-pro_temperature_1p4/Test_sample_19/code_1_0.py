def music_analysis():
    """
    This function provides a step-by-step analysis of the modulation
    in Beethoven's Moonlight Sonata as requested.
    """
    part1_answer = "1. To what key the music modulates in measure 11 and the beginning of measure 12?\n   The music modulates to the key of A minor. The G-sharp diminished chord in measure 11 functions as the leading-tone chord to A minor, resolving to the A minor triad on the first beat of measure 12."
    
    part2_answer = "2. Given the 4 sharps key environment, what is the connection and / or justification for this particular modulation?\n   The original key is C-sharp minor. The key of A is the submediant (VI) of C-sharp minor. Typically, the submediant chord would be A major. By modulating to A minor instead, Beethoven uses 'modal mixture' (borrowing from the parallel minor), which creates a more somber and harmonically distant effect. This modulation can also be viewed as moving to the subdominant (iv) of the relative major key (E major)."

    part3_answer = "3. What would be the complete and accurate Roman numeral marking for the first beat of measure 11?\n   The chord is a G-sharp diminished triad, which is the leading-tone chord (vii°) of the new key, A. Since A is the submediant (VI) of the home key (C-sharp minor), the chord is analyzed as a secondary or applied chord. The correct Roman numeral is vii°/VI."

    print(part1_answer)
    print("\n" + "-"*20 + "\n")
    print(part2_answer)
    print("\n" + "-"*20 + "\n")
    print(part3_answer)

music_analysis()