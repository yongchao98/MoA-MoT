def solve_task():
    """
    This function provides the solution by identifying the author and analyzing the text.
    """

    # 1. Identify the author.
    # The text is from Chapter 75 of Maimonides' "Guide for the Perplexed".
    # The author's name is Maimonides.
    author_name = "maimonides"

    # 2. Identify the stressed syllables for the first 10 words.
    # The first 10 words in Judaeo-Arabic are:
    # אנא אכין לך איצא פי הרא אלפצל ולאיל אלתוחיד עלי
    # Transcribed into Arabic:
    # أنا أبين لك أيضا في هذا الفصل دلائل التوحيد على
    # (anā ubayyinu laka ayḍan fī hādhā al-faṣli dalāʾila al-tawḥīdi ʿalā)
    #
    # Applying Modern Standard Arabic stress rules:
    # 1. أنا (ʔa-nā) -> stressed on the first syllable (exceptional) -> 'a (ا)
    # 2. أبين (ʔu-bay-yi-nu) -> stressed on the antepenult -> bay (بي)
    # 3. لك (la-ka) -> stressed on the first syllable -> la (ل)
    # 4. أيضا (ʔay-ḍan) -> stressed on the penult -> 'ay (اي)
    # 5. في (fī) -> monosyllabic -> fī (في)
    # 6. هذا (hā-dhā) -> stressed on the penult -> hā (ها)
    # 7. الفصل (al-faṣ-li) -> stress on the root word's penult -> faṣ (فص)
    # 8. دلائل (da-lā-ʔi-li) -> stressed on the antepenult -> lā (لا)
    # 9. التوحيد (at-taw-ḥīd) -> stressed on the superheavy final syllable -> ḥīd (حيد)
    # 10. على (ʿa-lā) -> stressed on the final syllable -> lā (لى)
    #
    # The list of stressed syllables in Arabic script without vocalization:
    stressed_syllables = "ا بي ل اي في ها فص لا حيد لى"

    # Print the final answer in the required format.
    print(author_name)
    print()  # Exactly one line space
    print(stressed_syllables)

solve_task()