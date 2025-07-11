def analyze_moonlight_sonata():
    """
    Analyzes a passage from Beethoven's Moonlight Sonata and prints the findings.
    """

    # --- Part 1: Key Identification ---
    # The chord in m. 11 (G-B-D) acts as a dominant (V) to the chord
    # in m. 12 (C-E-G), which acts as a tonic (I).
    modulated_key = "C Major"

    # --- Part 2: Justification for Modulation ---
    original_key = "C# minor"
    relative_major_key = "E Major"
    # The new key, C Major, is the flattened submediant (bVI) of the relative major.
    # This is a chromatic modulation common for the era.
    relationship_explanation = (
        f"The original key is {original_key}. Its relative major is {relative_major_key}.\n"
        f"The modulated key, {modulated_key}, is the flattened submediant (bVI) of the relative major ({relative_major_key}).\n"
        "This is a dramatic chromatic modulation, creating a striking color change."
    )

    # --- Part 3: Roman Numeral Analysis ---
    # The chord in m. 11 is G-B-D, which is G Major in first inversion.
    # It functions as a secondary dominant to the new key area.
    numeral_function = "V"
    inversion_figure = "6"
    # The destination key (C Major) is analyzed as bVI of the relative major (E Major).
    destination_numeral = "bVI"
    complete_numeral = f"{numeral_function}⁶/{destination_numeral}"

    # --- Print the full analysis ---
    print("--- Analysis of Moonlight Sonata, m. 11-12 ---")

    print("\nPart 1: To what key the music modulates?")
    print(f"The music modulates to the key of {modulated_key}.")

    print("\nPart 2: What is the connection and/or justification for this modulation?")
    print(relationship_explanation)

    print("\nPart 3: What is the complete and accurate Roman numeral marking for the first beat of measure 11?")
    print("The chord is the dominant of the new key area, in first inversion.")
    print(f"Its function is the dominant, represented by the Roman numeral: {numeral_function}")
    print(f"It is in first inversion, represented by the figure: {inversion_figure}")
    print(f"It tonicizes the destination key area, which is the flattened submediant, represented by: {destination_numeral}")
    print(f"The complete Roman numeral is: {complete_numeral}")


if __name__ == "__main__":
    analyze_moonlight_sonata()
    final_answer = "1. C Major. 2. The modulation is to the key of the flattened submediant (♭VI) of the piece's relative major key (E Major). This is a form of chromatic modulation. 3. V⁶/♭VI"
    print(f"\n<<< {final_answer} >>>")
