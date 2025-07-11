def get_roman_numeral():
    """
    This function analyzes the chord in measure 30 of Beethoven's "Pathetique"
    sonata (1st mvt.) and provides the corresponding Roman numeral.
    """
    key = "C minor"
    measure = 30
    
    # Notes on the first beat of measure 30
    notes = ["G", "B-natural", "F"]
    
    # Analysis
    root = "G"
    scale_degree = "5th"
    scale_degree_roman = "V"
    chord_quality = "Dominant 7th"
    extension = "7"
    
    print(f"Analysis of Beethoven's 'Pathetique' Sonata, 1st Mvt., Measure {measure}:")
    print(f"1. The key is {key}.")
    print(f"2. The notes in the chord are: {', '.join(notes)}.")
    print(f"3. The root of the chord is {root}, which is the {scale_degree} degree of the {key} scale.")
    print(f"4. The chord is a {chord_quality}, which is represented by a capital Roman numeral for the scale degree and the number for the seventh.")
    
    print("\nConstructing the final Roman numeral:")
    print(f"The scale degree is the Dominant: {scale_degree_roman}")
    print(f"The chord extension is the seventh: {extension}")
    
    final_numeral = scale_degree_roman + extension
    print(f"\nTherefore, the correct Roman numeral is {final_numeral}.")

get_roman_numeral()
<<<V7>>>