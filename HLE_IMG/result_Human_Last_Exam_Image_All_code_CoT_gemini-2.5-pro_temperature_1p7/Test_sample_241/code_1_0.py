def solve_maimonides_text():
  """
  Solves the task by identifying the author and analyzing word stress in the provided text.
  """
  # 1. Author Identification
  # The text is from Part 1, Chapter 75 of "The Guide for the Perplexed".
  author_name = "maimonides"

  # 2. Stressed Syllable Analysis
  # The first 10 words of the main body are transcribed into Arabic and analyzed for stress
  # based on Modern Standard Arabic rules.
  # 1. anā (أَنَا) -> 'a-nā -> Stress on 'a (special case) -> أ
  # 2. ubayyinu (أُبَيِّنُ) -> u-bay-yi-nu -> Stress on bay -> بي
  # 3. laka (لَكَ) -> la-ka -> Stress on la -> ل
  # 4. ayḍan (أَيْضًا) -> 'ay-ḍan -> Stress on 'ay -> اي
  # 5. fī (فِي) -> fī -> Stress on fī -> في
  # 6. hādhā (هَٰذَا) -> hā-dhā -> Stress on hā -> ها
  # 7. al-faṣli (الْفَصْلِ) -> 'al-faṣ-li -> Stress on faṣ -> فص
  # 8. wa-dalāʾila (وَدَلَائِلَ) -> wa-da-lā-'i-la -> Stress on lā -> لا
  # 9. al-tawḥīdi (التَّوْحِيدِ) -> 'at-taw-ḥī-di -> Stress on ḥī -> حي
  # 10. ʿalā (عَلَىٰ) -> 'a-lā -> Stress on lā -> لى
  stressed_syllables = ["أ", "بي", "ل", "اي", "في", "ها", "فص", "لا", "حي", "لى"]

  # Print the final output as specified.
  print(author_name)
  print() # Leave exactly one line space
  print(" ".join(stressed_syllables))

solve_maimonides_text()
#<<<maimonides
#
#أ بي ل اي في ها فص لا حي لى>>>