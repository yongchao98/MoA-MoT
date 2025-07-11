def analyze_harmony():
    """
    Analyzes the harmony of measure 30 of Beethoven's "Pathetique" Sonata, 1st movement.
    """
    
    # Step 1: Identify the notes in the measure
    notes_in_measure_30 = ["A", "C#", "E"]
    chord_in_measure_30 = "A Major"

    # Step 2: Identify the chord in the following measure
    notes_in_measure_31 = ["D", "F", "A"]
    chord_in_measure_31 = "D minor"

    # Step 3: Define the home key and the relationship between the chords
    home_key = "C minor"
    supertonic_chord = "D minor (ii)"
    relationship = "A Major is the dominant (V) of D minor."

    # Step 4: Determine the Roman Numeral
    roman_numeral = "V/ii"

    # Print the analysis
    print("Analysis of Beethoven 'Pathetique' Sonata, 1st mvt., measure 30:")
    print("-" * 60)
    print(f"1. The notes in measure 30 are {', '.join(notes_in_measure_30)}, which form an {chord_in_measure_30} chord.")
    print(f"2. The home key of the movement is {home_key}.")
    print(f"3. The chord in the next measure (m. 31) is {chord_in_measure_31}.")
    print(f"4. In the key of {home_key}, the {supertonic_chord} chord is the supertonic (ii).")
    print(f"5. The {chord_in_measure_30} chord functions as the dominant of the {supertonic_chord} chord.")
    print("6. This type of chord is called a 'secondary dominant'.")
    print("-" * 60)
    print(f"Therefore, the correct Roman numeral is written as 'V of ii'.")
    print(f"Final Answer: {roman_numeral}")

analyze_harmony()