def get_roman_numeral_analysis():
    """
    Analyzes the harmony of measure 30 of Beethoven's Pathetique Sonata, 1st mvt.
    and prints the corresponding Roman numeral.
    """

    key_signature = "C minor"
    measure_number = 30

    # Notes in measure 30
    left_hand_notes = ["G", "B-natural", "D"]
    right_hand_notes = ["F-sharp", "G", "A-sharp", "B-natural"]
    underlying_chord = "G Major"

    # Analysis
    root_of_chord = "G"
    key_tonic = "C"
    scale_degree = 5  # G is the 5th degree of the C minor scale

    # In Roman numeral analysis, a major chord on the dominant (5th degree) is 'V'
    roman_numeral = "V"

    print(f"Analysis of Beethoven's 'Pathetique' Sonata, 1st mvt., measure {measure_number}:")
    print(f"1. The key is {key_signature}.")
    print(f"2. The notes in the left hand are {', '.join(left_hand_notes)}, forming a {underlying_chord} chord.")
    print(f"3. In the key of {key_signature}, the note '{root_of_chord}' is the {scale_degree}th scale degree (the dominant).")
    print(f"4. A major chord built on the dominant is represented by the Roman numeral '{roman_numeral}'.")
    print(f"\nThe correct Roman numeral for measure {measure_number} is: {roman_numeral}")

get_roman_numeral_analysis()