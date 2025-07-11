def analyze_beethoven_measure():
    """
    Analyzes measure 30 of Beethoven's "Pathetique" Sonata, Mvt. 1,
    to determine its Roman numeral.
    """
    # Step 1: Define the musical context.
    local_key = "E-flat Major"
    target_chord_in_key = "C minor"
    target_chord_roman_numeral = "vi"

    # Step 2: Identify the chord in measure 30.
    chord_m30_notes = "G, B-natural, D, F"
    chord_m30_name = "G dominant 7th"
    chord_m30_root = "G"

    # Step 3: Determine the function.
    # The G7 chord is the dominant of the C minor chord.
    # Since C minor is the 'vi' in E-flat Major, G7 is the 'V of vi'.
    function_description = f"V of {target_chord_roman_numeral}"
    
    # Step 4: Assemble the final Roman numeral components.
    dominant_symbol = "V"
    seventh_symbol = "⁷"
    separator_symbol = "/"
    target_symbol = "vi"
    
    final_numeral = f"{dominant_symbol}{seventh_symbol}{separator_symbol}{target_symbol}"

    print("Analysis Steps:")
    print(f"1. The local key is {local_key}.")
    print(f"2. The chord in measure 30 is a {chord_m30_name} ({chord_m30_notes}).")
    print(f"3. This chord resolves to {target_chord_in_key}, which is the '{target_chord_roman_numeral}' chord in {local_key}.")
    print(f"4. Thus, the chord functions as a secondary dominant: the {function_description}.")
    
    print("\nFinal Roman Numeral Construction:")
    print(f"Dominant function: {dominant_symbol}")
    print(f"Seventh quality: {seventh_symbol}")
    print(f"Secondary function ('of'): {separator_symbol}")
    print(f"Target chord: {target_symbol}")

    print(f"\nThe complete Roman numeral is: {final_numeral}")


analyze_beethoven_measure()
<<<V⁷/vi>>>