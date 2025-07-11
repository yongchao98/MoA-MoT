# -*- coding: utf-8 -*-

def solve_hebrew_linguistics():
    """
    This function solves the two-part linguistic task and prints the result.
    Part 1: Identifies primary stressed syllables in the first seven words of Psalm 74:1.
    Part 2: Identifies the last syllable with marked secondary stress in 1 Chronicles 5:10.
    The final output string combines both answers as requested.
    """
    
    # Analysis for Part 1: Primary stressed syllables
    # 1. מַשְׂכִּ֗יל -> Stress on כִּיל (marked by Rivia) -> כִּ֗יל
    # 2. לְאָ֫סָ֥ף -> Stress on אָ (marked by Ole) -> אָ֫
    # 3. לָמָ֣ה -> Stress on מָה (marked by Munach) -> מָ֣ה
    # 4. אֱ֭לֹהִים -> Stress on לֹ (marked by Tipeha) -> לֹ֭
    # 5. זָנַ֣חְתָּ -> Stress on נַח (marked by Munach) -> נַ֣ח
    # 6. לָנֶ֑צַח -> Stress on נֶ (marked by Etnachta) -> נֶ֑
    # 7. יֶעְשַׁ֥ן -> Stress on שַׁן (marked by Merkha) -> שַׁ֥ן
    answer_part_1 = "כִּ֗יל אָ֫ מָ֣ה לֹ֭ נַ֣ח נֶ֑ שַׁ֥ן"
    
    # Analysis for Part 2: Last marked secondary stress syllable
    # Verse: וּבִימֵ֣י שָׁא֗וּל עָשׂ֤וּ מִלְחָמָה֙ עִם־הַֽהַגְרִאִ֔ים וַֽיִּפְּל֖וּ בְּיָדָ֑ם וַיֵּשְׁבוּ֙ בְּאָ֣הֳלֵיהֶ֔ם עַֽל־כָּל־פְּנֵ֖י מִזְרָ֥ח לַגִּלְעָֽד׃
    # The last meteg (marking secondary stress) appears on the word עַל־.
    # The syllable is עַל -> עַֽל
    answer_part_2 = "עַֽל"
    
    # Combine the answers as per the formatting rules.
    final_answer = f"{answer_part_1},{answer_part_2}"
    
    print(final_answer)

solve_hebrew_linguistics()