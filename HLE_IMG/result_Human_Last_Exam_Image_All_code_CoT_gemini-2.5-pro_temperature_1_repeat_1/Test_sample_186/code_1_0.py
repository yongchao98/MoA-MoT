def solve_manuscript_task():
    """
    This function provides the solution to the manuscript analysis task.
    
    Task 1: Identify the Verse.
    The text 'وياسوبو موسا وهرون إل فرعون ويومر لهم' corresponds to the Hebrew
    'וַיָּשֻׁבוּ מֹשֶׁה וְאַהֲרֹן אֶל-פַּרְעֹה וַיֹּאמֶר לָהֶם', which is a textual variant of Exodus 10:8.
    The answer is Exo. 10:8.
    
    Task 2: Compare Matres Lectionis.
    The comparison is done word by word, right to left.
    - 'وياسوبو' (וַיָּשֻׁבוּ): Arabic adds 'ا' for the 'a' vowel in 'יָ'. Result: ا
    - 'موسا' (מֹשֶׁה): Arabic adds 'و' for the 'o' vowel. Hebrew 'ה' is substituted by Arabic 'ا'. Result: و הا
    - 'فرعون' (פַּרְעֹה): Hebrew 'ה' is substituted by Arabic 'و'. Result: הو
    - 'ويومر' (וַיֹּאמֶר): Hebrew 'א' is substituted by Arabic 'و'. Result: או
    
    The final combined answer is formatted as requested.
    """
    
    # The verse identification for Task 1
    verse_id = "Exo. 10:8"
    
    # The comparison of matres lectionis for Task 2, in order of appearance (right to left)
    comparison_results = "ا و הا הو או"
    
    # Combine the answers with a single comma and no space
    final_answer = f"{verse_id},{comparison_results}"
    
    print(final_answer)

solve_manuscript_task()