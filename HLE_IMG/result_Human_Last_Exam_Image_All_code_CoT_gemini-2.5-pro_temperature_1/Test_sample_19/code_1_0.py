def solve_music_theory():
    """
    This function provides the answers to the three-part music theory question
    about Beethoven's Moonlight Sonata.
    """
    answer_part_1 = "1. Key of Modulation: The music modulates to the key of A minor."
    
    answer_part_2 = ("2. Justification: The original key is C# minor. Its relative major is E major. "
                     "The new key, A minor, is the subdominant minor (iv) of E major. "
                     "This modulation to the subdominant of the relative major is a common and effective harmonic device.")
                     
    answer_part_3 = ("3. Roman Numeral: The complete and accurate Roman numeral marking for the "
                     "first beat of measure 11 is vii°/vi. This represents the secondary "
                     "leading-tone diminished triad of the minor submediant (A minor in the key of C# minor).")

    full_answer = f"{answer_part_1}\n{answer_part_2}\n{answer_part_3}"
    
    final_answer_formatted = "<<<1. A minor\n2. The new key, A minor, is the subdominant minor (iv) of the original key's relative major (E Major).\n3. vii°/vi>>>"

    print(full_answer)
    print(f"\nFinal Answer Summary:\n{final_answer_formatted}")

solve_music_theory()