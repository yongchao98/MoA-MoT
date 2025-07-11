def analyze_measure_30():
    """
    Analyzes the harmony of measure 30 from Beethoven's Pathetique Sonata, 1st movement,
    and provides the correct Roman numeral.
    """
    key = "C minor"
    measure = 30
    notes_in_chord = ["G-sharp", "B-natural", "D", "F"]
    chord_name = "G-sharp fully diminished seventh (G#°7)"

    print(f"Analysis of Beethoven's Pathetique Sonata, 1st Movement, Measure {measure}:")
    print(f"1. Home Key: {key}")
    print(f"2. Notes in measure {measure}: {', '.join(notes_in_chord)}.")
    print(f"3. Identified Chord: {chord_name}.")
    print("\nFunctional Analysis:")
    print("   - The root of the chord, G-sharp, is the leading tone to the note A.")
    print("   - This chord functions as a secondary leading-tone chord, tonicizing the submediant (the 6th degree).")
    print("   - In C minor, the diatonic submediant is VI (A-flat major), but Beethoven is targeting the borrowed minor submediant, vi (A minor).")
    print("\nResulting Roman Numeral Breakdown:")
    print("   - 'vii' represents a leading-tone triad.")
    print("   - The superscript '°7' indicates it is a fully diminished seventh chord.")
    print("   - The '/vi' indicates it functions as the leading-tone chord 'of the minor submediant' (the 6th degree chord).")
    
    # As requested, outputting the numbers in the final notation
    numeral_7 = 7
    numeral_vi = "vi" # In Roman numerals, 6 is VI or vi
    numeral_6_from_vi = 6
    
    print("\nThe final Roman numeral is vii°7/vi.")
    print(f"The first number in the notation is {numeral_7}.")
    print(f"The number represented by 'vi' is {numeral_6_from_vi}.")

analyze_measure_30()