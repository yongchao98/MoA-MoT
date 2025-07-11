def solve_mozart_analysis():
    """
    This function provides the Roman numeral analysis for the specified chord.
    """
    # Components of the Roman numeral
    root_function = "vii"
    inversion_figure_6 = "\u2076"  # Superscript 6
    inversion_figure_4 = "\u2074"  # Superscript 4
    secondary_function_target = "/iv"

    # Assemble the final Roman numeral
    roman_numeral = root_function + inversion_figure_6 + inversion_figure_4 + secondary_function_target

    print("The key is D minor.")
    print("The notes in the circled chord are C#, F#, and A, which form an F# minor triad.")
    print("The bass note is C#, the fifth of the chord, so it is in second inversion (⁶⁴).")
    print("The chord resolves to G minor, the subdominant (iv) chord.")
    print("Therefore, its function is the secondary leading-tone chord of the subdominant.")
    print("\nThe final Roman numeral is:")
    print(roman_numeral)

solve_mozart_analysis()