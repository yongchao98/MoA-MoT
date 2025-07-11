def find_roman_numeral():
    """
    Analyzes the harmony of measure 30 from Beethoven's "Pathetique" Sonata
    and prints the corresponding Roman numeral.
    """
    key = "C minor"
    measure = 30
    
    # Notes identified in the measure
    bass_note = "G"
    chord_notes = ["G", "B-natural", "D"]
    non_chord_tone = "E-flat"
    
    # Harmonic Analysis
    # The G major chord (G-B-natural-D) is the dominant of C minor.
    # The dominant chord is represented by the Roman numeral V.
    # It is an uppercase V because the chord is major.
    roman_numeral = "V"
    
    print(f"Step-by-step analysis for measure {measure}:")
    print(f"1. The key of the piece is {key}.")
    print(f"2. The bass note is {bass_note}, which is the dominant of {key}.")
    print(f"3. The essential chord tones are {', '.join(chord_notes)}, forming a G major triad.")
    print(f"4. This G major triad functions as the dominant chord in the key of {key}.")
    print("5. The Roman numeral for a major dominant chord is an uppercase V.")
    
    print("\nFinal Roman Numeral:")
    # The prompt asks to output each 'number' in the final 'equation'.
    # In this context, the Roman numeral 'V' is the final result.
    print(roman_numeral)

find_roman_numeral()