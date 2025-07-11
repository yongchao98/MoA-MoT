def analyze_mozart_chord():
    """
    Analyzes the chord in measure 8 of Mozart's Fantasy K. 397
    and provides its Roman numeral.
    """
    
    # Define the components of the analysis
    key = "D minor"
    notes_in_chord = ["F#", "A", "C#", "E"]
    target_chord_numeral = "iv"
    
    # Step 1: Determine the function based on the bass note
    functional_root = "F#"
    target_root = "G"
    function_numeral = "vii"
    
    # Step 2: Determine the chord quality
    # F#-A (m3), F#-C# (P5), F#-E (m7) -> Minor Seventh Chord
    quality = "m7"
    
    # Step 3: Combine the components for the full numeral
    full_numeral = f"{function_numeral}{quality}/{target_chord_numeral}"
    
    # Print the explanation and the final answer
    print("This program determines the Roman numeral for the chord in measure 8.")
    print(f"The key is {key}, and the chord leads to the subdominant ({target_chord_numeral}).")
    print("-" * 20)
    print("Analysis Breakdown:")
    print(f"1. The chord functions as a secondary leading-tone chord to 'iv'. Its numeral root is: {function_numeral}")
    print(f"2. The quality of the chord ({', '.join(notes_in_chord)}) is a minor seventh. Its quality symbol is: {quality}")
    print(f"3. The target of this secondary chord is the subdominant. This is represented as: /{target_chord_numeral}")
    print("-" * 20)
    print("The final Roman numeral is constructed by combining these parts:")
    print(f"Final Accurate Roman Numeral: {full_numeral}")

analyze_mozart_chord()