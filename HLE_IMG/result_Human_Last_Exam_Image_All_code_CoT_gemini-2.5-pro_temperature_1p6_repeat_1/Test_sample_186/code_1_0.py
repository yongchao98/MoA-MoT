def solve_manuscript_task():
    """
    This function provides the solution to the manuscript analysis task.
    1. Identifies the biblical verse from the specified lines.
    2. Compares the matres lectionis between the Hebrew BHS text and the Arabic transcription.
    3. Formats and prints the final answer as a single string.
    """
    
    # Part 1: Identify the Verse
    # The text corresponds to the Hebrew:
    # וְאֵלֶּה שְׁמוֹת בְּנֵי יִשְׂרָאֵל הַבָּאִים מִצְרָיְמָה אֵת יַעֲקֹב אִישׁ וּבֵיתוֹ בָּאוּ
    # This is Exodus 1:1.
    verse_identification = "Ex. 1:1"

    # Part 2: Compare Matres Lectionis
    # The analysis results are compiled in reading order (right-to-left).
    # וְאֵלֶּה -> والي : additional ا, substitution הي
    # שְׁמוֹת -> اسموت : match
    # בְּנֵי -> بني : match
    # יִשְׂרָאֵל -> يسرايل : additional ا, additional ي
    # הַבָּאִים -> البايين : additional ا
    # מִצְרָיְמָה -> مصريما : substitution הا
    # אֵת -> ايت : additional ي
    # יַעֲקֹב -> يعقوب : additional و
    # אִישׁ -> ايش : match
    # וּבֵיתוֹ -> وبيتو : match
    # בָּאוּ -> باوا : additional ا, substitution وا
    
    matres_comparison = [
        "ا הي",  # for וְאֵלֶּה
        "ا ي",   # for יִשְׂרָאֵל
        "ا",     # for הַבָּאִים
        "הا",    # for מִצְרָיְמָה
        "ي",     # for אֵת
        "و",     # for יַעֲקֹב
        "ا וا"   # for בָּאוּ
    ]
    
    # The final string is the concatenation of these findings, separated by spaces.
    matres_string = " ".join(matres_comparison)

    # Combine Part 1 and Part 2 with a single comma.
    final_answer = f"{verse_identification},{matres_string}"
    
    print(final_answer)

solve_manuscript_task()