def solve_music_theory_analysis():
    """
    This function analyzes the provided musical excerpt from Beethoven's Moonlight Sonata
    and prints the answers to the three-part question.
    """

    # --- Part 1: Identify the key of the modulation ---
    part1_answer = "1. To what key the music modulates in measure 11 and the beginning of measure 12?"
    part1_solution = ("The music modulates to the key of B major. The chord in measure 11 is E major, "
                      "which acts as the subdominant (IV) to the B major chord, the new tonic (I), "
                      "at the beginning of measure 12.")

    # --- Part 2: Justify the modulation ---
    part2_answer = ("2. Given the 4 sharps key environment, what is the connection and / or justification for this "
                    "particular modulation?")
    part2_solution = ("The home key of the piece is C# minor. Its relative major is E major. "
                      "The new key, B major, is the dominant (V) of E major. Therefore, the justification is that "
                      "the music moves from the tonic minor to the dominant of the relative major. "
                      "In the key of C# minor, this is a progression from i to V/III, which is a common and logical modulation.")

    # --- Part 3: Determine the Roman numeral marking ---
    part3_answer = ("3. Building on 1 and 2 answers - what would be the complete and accurate roman numeral marking "
                    "for the first beat of measure 11?")
    part3_solution = ("The complete and accurate Roman numeral marking for the first beat of measure 11 is IV/bVII. "
                      "Here's the breakdown:\n"
                      "   - The chord itself is E major.\n"
                      "   - Its function is the subdominant (IV) of the new key, B major.\n"
                      "   - The new key of B major is the bVII (built on the lowered 7th scale degree, or subtonic) "
                      "relative to the home key of C# minor.\n"
                      "   - Thus, 'IV/bVII' means 'the subdominant of the key of the lowered subtonic'.")
    
    # --- Print the final results ---
    print(f"{part1_answer}\nANSWER: {part1_solution}\n")
    print(f"{part2_answer}\nANSWER: {part2_solution}\n")
    print(f"{part3_answer}\nANSWER: {part3_solution}")

solve_music_theory_analysis()