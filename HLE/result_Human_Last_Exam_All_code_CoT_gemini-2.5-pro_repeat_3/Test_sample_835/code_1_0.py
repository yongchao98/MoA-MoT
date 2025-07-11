def solve_hebrew_tasks():
    """
    This function solves the two-part Hebrew analysis task and prints the result.
    """
    # Part 1: List of syllables with primary word stress.
    primary_stress_syllables = "כִּיל אָ לָ לֹ נַח לָ שַׁן"
    
    # Part 2: Last occurrence of a syllable with marked secondary stress.
    last_secondary_stress_syllable = "לַ"
    
    # Combine the results as per the specified format.
    # Format: [Part 1 result],[Part 2 result] with no space after the comma.
    final_answer = f"{primary_stress_syllables},{last_secondary_stress_syllable}"
    
    print(final_answer)

solve_hebrew_tasks()