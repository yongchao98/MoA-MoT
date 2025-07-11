def solve_music_theory():
    """
    This function provides the answers to the three-part music theory question
    about Beethoven's Moonlight Sonata.
    """
    
    # Answer to Part 1
    key_modulation = "A minor"
    
    # Answer to Part 2
    justification = (
        "The home key is C-sharp minor. The new key, A minor, is the minor submediant (vi) of C-sharp minor. "
        "The submediant (the 6th scale degree) is a closely related key area, making this a smooth and logical, "
        "yet colorful, chromatic modulation."
    )
    
    # Answer to Part 3
    roman_numeral = "vii°7/vi"
    
    print("--- Moonlight Sonata Analysis ---")
    print("\nPart 1: To what key the music modulates in measure 11 and the beginning of measure 12?")
    print(f"Answer: The music modulates to the key of {key_modulation}.")
    
    print("\nPart 2: Given the 4 sharps key environment, what is the connection and / or justification for this particular modulation?")
    print(f"Answer: {justification}")
    
    print("\nPart 3: What would be the complete and accurate roman numeral marking for the first beat of measure 11?")
    print(f"Answer: The Roman numeral is {roman_numeral}. This represents the leading-tone seventh chord of the minor submediant.")

solve_music_theory()
# The final combined answer is constructed from the analysis above.
final_answer = (
    "1. The music modulates to A minor. "
    "2. The home key is C# minor. A is the submediant (6th degree). The modulation is to the minor submediant (vi), a common chromatic modulation to a closely related key. "
    "3. The Roman numeral is vii°7/vi."
)
print(f"\n<<<1. The music modulates to A minor. 2. The home key is C# minor. A is the submediant (6th degree). The modulation is to the minor submediant (vi), a common chromatic modulation to a closely related key. 3. The Roman numeral is vii°7/vi.>>>")