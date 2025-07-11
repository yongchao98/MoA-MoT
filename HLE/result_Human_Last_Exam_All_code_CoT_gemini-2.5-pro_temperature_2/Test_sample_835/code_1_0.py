def solve_hebrew_tasks():
    """
    This function generates the solution to the two Hebrew linguistics tasks.
    """
    # Solution for task 1: Primarily stressed syllables in the first seven words.
    # 1. מַשְׂכִּ֗יל -> כִּ֗יל
    # 2. לְאָ֫סָף -> אָ֫
    # 3. לָמָ֣ה -> לָ֣
    # 4. אֱ֭לֹהִים -> לֹ֭ (syllable with the te'amim mark)
    # 5. זָנַ֣חְתָּ -> נַ֣ח
    # 6. לָנֶ֑צַח -> נֶ֑
    # 7. יֶעְשַׁ֥ן -> שַׁ֥ן
    stressed_syllables = "כִּ֗יל אָ֫ לָ֣ לֹ֭ נַ֣ח נֶ֑ שַׁ֥ן"

    # Solution for task 2: Last occurrence of a marked secondary stress.
    # Scanning from the end of the verse, the last word is לַגִּלְעָֽד׃.
    # This word contains a meteg on the syllable עָד.
    secondary_stress_syllable = "עָֽד"

    # Combine the results as per the user's formatting instructions.
    final_answer = f"{stressed_syllables},{secondary_stress_syllable}"
    
    print(final_answer)

solve_hebrew_tasks()