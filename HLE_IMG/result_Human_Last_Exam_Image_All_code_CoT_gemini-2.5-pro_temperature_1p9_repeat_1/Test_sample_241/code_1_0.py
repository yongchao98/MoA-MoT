def solve_task():
    """
    This script provides the solution by identifying the author of the text
    and the primary stressed syllables of the first 10 words.
    """
    # 1. The author of the book is Moses Maimonides.
    author = "maimonides"

    # 2. The first 10 words are:
    # Judaeo-Arabic: אנא אכין לך איצא פי הרא אלפצל ולאיל אלתוחיד עלי
    # Arabic: أنا أبين لك أيضا في هذا الفصل ودلائل التوحيد على
    # The list of syllables with primary word stress is derived from applying
    # Modern Standard Arabic stress rules to these words.
    # anā -> 'nā' -> نا
    # ubayyinu -> 'bay' -> بي
    # laka -> 'la' -> ل
    # ayḍan -> 'ay' -> اي
    # fī -> 'fī' -> في
    # hādhā -> 'hā' -> ها
    # al-faṣlu -> 'faṣ' -> فص
    # wa-dalāʾilu -> 'lā' -> لا
    # al-tawḥīdi -> 'ḥī' -> حي
    # ʿalā -> 'lā' -> لا
    stressed_syllables = "نا بي ل اي في ها فص لا حي لا"
    
    # Printing the result in the required format.
    print(author)
    print()
    print(stressed_syllables)

solve_task()