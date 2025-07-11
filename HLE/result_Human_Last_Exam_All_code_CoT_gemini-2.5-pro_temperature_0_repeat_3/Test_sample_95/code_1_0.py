def solve_riddle():
    """
    This script breaks down the riddle's clues to find the answer.
    """
    # The clues from the riddle
    clue_location = "Smoky northern cities vs. Milan"
    clue_time = "19th century"
    clue_subject = "Something impossible to see"
    clue_person = "Kasimir Graf"
    clue_reason = "Lack of imagination"

    print("--- Riddle Analysis ---")
    print(f"1. Location Contrast ('{clue_location}'):")
    print("   - In the 19th century, industrial cities in northern Europe had smoggy skies, terrible for sky-watching.")
    print("   - Milan, Italy, had clearer atmospheric conditions, ideal for astronomy.")
    print("\n" + "=" * 25 + "\n")

    print(f"2. Historical Context ('{clue_time}'):")
    print("   - This was the peak of public interest in Mars, sparked by Italian astronomer Giovanni Schiaparelli.")
    print("   - From Milan, he observed 'canali' (channels) on Mars, which was mistranslated as 'canals'.")
    print("\n" + "=" * 25 + "\n")

    print(f"3. The Conclusion ('{clue_subject}', '{clue_person}', '{clue_reason}'):")
    print("   - The 'canals' led to the belief in intelligent life on Mars.")
    print("   - Seeing these faint features was difficult and controversial; many astronomers couldn't see them and claimed it required 'imagination'.")
    print("   - Therefore, in smoky cities, it was impossible to see the evidence for THEM.")
    print("\n" + "=" * 25 + "\n")

    # The final "equation" as requested
    print("Final Equation:")
    print("   Clear Skies in Milan")
    print("+  Schiaparelli's 'Canali' Observations")
    print("+  Widespread 'Canal' Speculation")
    print("-  Poor Visibility in Smoky Cities")
    print("---------------------------------------")
    print("=  THEM, who could only be seen with 'imagination' by some.")

    answer = "Martians"
    print(f"\nTherefore, THEM are: {answer}")

solve_riddle()