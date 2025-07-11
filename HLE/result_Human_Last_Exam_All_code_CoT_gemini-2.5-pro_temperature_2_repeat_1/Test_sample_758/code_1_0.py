def find_common_elements():
    """
    This function analyzes the common literary elements between
    "The Yellow Wallpaper" and "Key West" and prints the result.
    """
    
    # Roman numerals corresponding to the statements
    # I, II, III, IV, V, VI
    numerals = ["I", "II", "III", "IV", "V", "VI"]

    # Analysis for "The Yellow Wallpaper"
    # [Confinement, Detachment, Indifference, Reunion Attempt, External Control, Medical]
    yellow_wallpaper_truths = [True, True, True, True, True, True]

    # Analysis for "Key West"
    # [Confinement, Detachment, Indifference, Reunion Attempt, External Control, Medical]
    key_west_truths = [False, True, True, True, True, False]

    common_features = []
    # Identify which numerals correspond to elements true for both stories
    for i in range(len(numerals)):
        if yellow_wallpaper_truths[i] and key_west_truths[i]:
            common_features.append(numerals[i])

    # Format the final list as a string, outputting each "number"
    # as requested in the prompt.
    final_answer = ", ".join(common_features)
    
    print(final_answer)

find_common_elements()