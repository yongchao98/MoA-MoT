def get_roman_numeral():
    """
    Analyzes and prints the Roman numeral for the chord in measure 30 of
    Beethoven's "Pathetique" Sonata, 1st movement.
    """
    key = "C minor"
    measure = 30
    
    # Step 1: Identify notes in the measure
    notes_in_chord = ["G", "B-natural", "D", "F-sharp"]
    
    # Step 2: Determine the chord's root and quality
    root = "G"
    quality = "Major 7th"
    chord_name = f"{root}{quality.replace(' ', '')}" # Gmaj7

    # Step 3: Determine the Roman numeral based on the key
    # G is the 5th scale degree in C minor, which is the dominant.
    roman_numeral_base = "V"
    
    # The chord has a major third (B-natural) and a major seventh (F-sharp)
    roman_numeral_quality = "maj7"
    final_roman_numeral = roman_numeral_base + roman_numeral_quality

    # Print the analysis step-by-step
    print(f"Analysis for measure {measure} of Beethoven's 'Pathetique' Sonata (1st Mov.):")
    print("-" * 70)
    print(f"1. The key is {key}.")
    print(f"2. The notes played in measure {measure} are: {', '.join(notes_in_chord)}.")
    print(f"3. These notes form a {root} {quality} chord ({chord_name}).")
    print(f"4. The root of the chord, {root}, is the dominant (5th degree) of {key}.")
    
    print("\nDeriving the Roman Numeral 'Equation':")
    print(f"  - The Roman numeral for the dominant is: {roman_numeral_base}")
    print(f"  - The quality of the chord (major triad with a major seventh) is represented by: {roman_numeral_quality}")
    
    print("\nFinal Result:")
    print(f"The correct Roman numeral is formed by combining these two parts.")
    print(f"Therefore, the final Roman numeral is: {final_roman_numeral}")

get_roman_numeral()