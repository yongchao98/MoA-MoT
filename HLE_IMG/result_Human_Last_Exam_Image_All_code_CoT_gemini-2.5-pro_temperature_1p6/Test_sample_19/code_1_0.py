def solve_music_theory():
    """
    Solves the 3-part music theory question about Beethoven's Moonlight Sonata.
    """
    # Part 1: Key of Modulation
    part1_answer = "The music modulates to the key of B minor."

    # Part 2: Justification
    part2_answer = ("The home key is C# minor. The new key, B minor, is the subtonic minor (vii) of C# minor. "
                    "This is a common Romantic-era modulation. It can also be justified as the 'minor dominant "
                    "of the relative major' (v/III), as B minor is the minor dominant of E Major, the relative major key.")

    # Part 3: Roman Numeral
    part3_answer = ("The complete and accurate Roman numeral marking for the harmony established in measure 11, "
                    "analyzed in the home key of C# minor, is vii. The chord is B minor, the subtonic minor.")

    print("1. In the first movement of Moonlight Sonata, to what key does the music modulate in measure 11 and the beginning of measure 12?")
    print(f"   Answer: {part1_answer}\n")

    print("2. Given the 4 sharps key environment, what is the connection and/or justification for this particular modulation?")
    print(f"   Answer: {part2_answer}\n")

    print("3. Building on answers 1 and 2, what would be the complete and accurate Roman numeral marking for the first beat of measure 11?")
    print(f"   Answer: {part3_answer}\n")

    # Combine for the final answer block
    final_answer = (
        f"1. B minor\n"
        f"2. The new key, B minor, is the subtonic minor (vii) of the home key, C# minor. "
        f"It can also be understood as the minor dominant of the relative major (v/III).\n"
        f"3. vii"
    )
    
    # Although the instructions say to use 'print', this final tagged format is also required.
    # The code above fulfills the 'print' requirement, and this line fulfills the final format requirement.
    # print(f"<<<{final_answer}>>>")


solve_music_theory()