def solve_hebrew_tasks():
    """
    This function solves the two Hebrew linguistics tasks and prints the result in the specified format.
    """

    # Task 1: Identify primary stressed syllables in the first seven words of the verse.
    # Verse: מַשְׂכִּ֗יל לְאָ֫סָ֥ף לָמָ֣ה אֱ֭לֹהִים זָנַ֣חְתָּ לָנֶ֑צַח יֶעְשַׁ֥ן אַ֝פְּךָ֗ בְּצֹ֣אן מַרְעִיתֶֽךָ׃
    # 1. מַשְׂכִּ֗יל -> Stress on כִּיל (marked by revia)
    # 2. לְאָ֫סָ֥ף -> Stress on אָ (marked by ole)
    # 3. לָמָ֣ה -> Stress on מָה (marked by mapakh)
    # 4. אֱ֭לֹהִים -> Ultimate stress on הִים (meteg marks secondary stress)
    # 5. זָנַ֣חְתָּ -> Stress on נַח (marked by zaqef qatan)
    # 6. לָנֶ֑צַח -> Stress on נֶ (marked by etnachta)
    # 7. יֶעְשַׁ֥ן -> Stress on שַׁן (marked by pashta)
    answer_1 = "כִּיל אָ מָה הִים נַח נֶ שַׁן"

    # Task 2: Identify the last occurrence of a syllable with secondary stress.
    # Verse: וּבִימֵ֣י שָׁא֗וּל עָשׂ֤וּ מִלְחָמָה֙ עִם־הַֽהַגְרִאִ֔ים וַֽיִּפְּל֖וּ בְּיָדָ֑ם וַיֵּשְׁבוּ֙ בְּאָ֣הֳלֵיהֶ֔ם עַֽל־כָּל־פְּנֵ֖י מִזְרָ֥ח לַגִּלְעָֽד׃
    # Secondary stress is marked by a meteg (ֽ). Scanning from the end of the verse,
    # the last word with a meteg is עַֽל־. The syllable is עַֽ.
    # The maqqef is to be omitted as per instructions.
    answer_2 = "עַֽ"

    # Combine the answers as per the required format: "Answer1,Answer2"
    final_answer = f"{answer_1},{answer_2}"

    # Print the final result.
    print(final_answer)

solve_hebrew_tasks()