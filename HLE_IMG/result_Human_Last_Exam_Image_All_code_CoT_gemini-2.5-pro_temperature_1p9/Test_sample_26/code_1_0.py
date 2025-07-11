def get_roman_numeral():
    """
    Analyzes the chord in measure 30 of Beethoven's "Pathetique" sonata,
    1st movement, and prints its Roman numeral analysis.
    """
    key = "C minor"
    bass_note = "G"
    right_hand_notes = ["D", "F#", "A", "C"]
    chord_name = "D7"
    function = "V7/V"

    print("--- Analysis of the Chord ---")
    print(f"The key of the piece is {key}.")
    print(f"The notes in measure 30 are {', '.join(right_hand_notes)} over a {bass_note} in the bass.")
    print(f"These notes form a {chord_name} chord.")
    print(f"In the key of {key}, the V (dominant) chord is G major.")
    print(f"The {chord_name} chord is the dominant of G, making it a secondary dominant.")
    print("Therefore, the function of this chord is the 'dominant of the dominant', or V of V.")

    print("\n--- Final Roman Numeral Equation ---")
    numeral_components = ['V', '7', '/', 'V']
    print(f"The final Roman numeral is: {''.join(numeral_components)}")
    print("\nPrinting each character of the numeral:")
    for char in numeral_components:
        print(char)

get_roman_numeral()