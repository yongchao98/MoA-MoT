def get_roman_numeral_analysis():
    """
    Provides a step-by-step harmonic analysis of measure 30
    of Beethoven's "Pathetique" Sonata, 1st movement.
    """

    key = "C minor"
    measure_chord_notes = "G#, C#, E"
    next_measure_chord = "D Major"

    print("Step 1: Identify the home key and the relevant chords.")
    print(f"The piece is in {key}.")
    print(f"The chord in measure 30 contains the notes: {measure_chord_notes}.")
    print(f"The chord in the next measure (m. 31) is {next_measure_chord}.\n")

    print("Step 2: Analyze the function of the chord in measure 31.")
    print(f"In the key of {key}, the dominant chord (V) is G Major.")
    print(f"The {next_measure_chord} chord is the dominant of G Major.")
    print("Therefore, the D Major chord is a secondary dominant, specifically the 'dominant of the dominant'.")
    print("Its Roman numeral is V/V.\n")

    print("Step 3: Analyze the function of the chord in measure 30.")
    print("The chord in measure 30 resolves to the D Major chord (V/V). Its function is to prepare the V/V.")
    print("The primary chords that prepare a harmony are its dominant or its leading-tone chord.")
    print("The leading tone of D is C#. The leading-tone triad of D is C# diminished (notes: C#-E-G).\n")

    print("Step 4: Compare the chord with its theoretical function.")
    print(f"The chord in measure 30 is {measure_chord_notes}.")
    print("This strongly resembles a C# diminished triad, with the G raised to G# for chromatic effect.")
    print("The G# can be understood as a chromatic passing tone from the G in measure 29.\n")

    print("Step 5: Determine the final Roman Numeral.")
    print("The chord's function is that of a leading-tone chord (vii°) applied to the V/V chord.")
    numeral_base = "vii°"
    applied_to = "V/V"
    final_numeral = f"{numeral_base}/{applied_to}"
    print(f"The Roman numeral is written as a leading-tone diminished triad '{numeral_base}' applied to the dominant of the dominant '{applied_to}'.")
    print("\nThe correct Roman numeral is:")
    print(f"'{numeral_base[0]}{numeral_base[1]}{numeral_base[2]}/{applied_to[0]}/{applied_to[2]}'")


get_roman_numeral_analysis()