def solve_task():
    """
    This function provides the solution by printing the author's name and the stressed syllables
    of the first 10 words of the provided text.
    """

    # 1. Author's name
    # The text is from "The Guide for the Perplexed" by Maimonides.
    author = "maimonides"

    # 2. List of stressed syllables for the first 10 words.
    # The words are: أنا أبين لك أيضا في هذا الفصل دلائل التوحيد على
    # Transliteration: anā ubayyinu laka ayḍan fī hādhā al-faṣli dalā'ila at-tawḥīdi ʿalā
    # Stress analysis:
    # 1. anā -> a-NĀ -> stressed syllable is 'nā' -> نا
    # 2. ubayyinu -> u-BAY-yi-nu -> stressed syllable is 'bay' -> بي
    # 3. laka -> LA-ka -> stressed syllable is 'la' -> ل
    # 4. ayḍan -> 'ay-ḌAN -> stressed syllable is 'ḍan' -> ضن
    # 5. fī -> FĪ -> stressed syllable is 'fī' -> في
    # 6. hādhā -> hā-DHĀ -> stressed syllable is 'dhā' -> ذا
    # 7. al-faṣli -> FAṢ-li -> stressed syllable is 'faṣ' -> فص
    # 8. dalā'ila -> da-LĀ-'i-la -> stressed syllable is 'lā' -> لا
    # 9. at-tawḥīdi -> taw-ḤĪ-di -> stressed syllable is 'ḥī' -> حي
    # 10. ʿalā -> ʿa-LĀ -> stressed syllable is 'lā' -> لا
    stressed_syllables = "نا بي ل ضن في ذا فص لا حي لا"

    # Print the final answer in the required format.
    print(author)
    print()  # For the required single blank line
    print(stressed_syllables)

solve_task()