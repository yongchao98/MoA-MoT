import sys
import io

# Ensure the output is printed in UTF-8 to handle Arabic characters.
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

def get_solution():
  """
  This function contains the solution to the problem.
  It identifies the author and determines the stressed syllables of the first ten words.
  """

  # Part 1: The name of the author of the book.
  # The text is from Part 1, Chapter 75 of "The Guide for the Perplexed", written by Maimonides.
  author_name = "maimonides"

  # Part 2: List of the syllables with primary word stress for the first 10 words.
  # The first 10 words are:
  # 1. ʾanā (أَنَا) -> a-NĀ -> Stress on 'nā'
  # 2. ʾubayyinu (أُبَيِّنُ) -> u-BAY-yi-nu -> Stress on 'bay'
  # 3. laka (لَكَ) -> LA-ka -> Stress on 'la'
  # 4. ʾayḍan (أَيْضاً) -> ay-ḌAN -> Stress on 'ḍan'
  # 5. fī (فِي) -> FĪ -> Stress on 'fī'
  # 6. hādhā (هَٰذَا) -> hā-DHĀ -> Stress on 'dhā'
  # 7. al-faṣli (الْفَصْلِ) -> al-FAṢ-li -> Stress on 'faṣ'
  # 8. wa-dalāʾila (وَدَلَائِلَ) -> wa-da-LĀ-ʾi-la -> Stress on 'lā'
  # 9. al-tawḥīdi (التَّوْحِيدِ) -> at-taw-ḤĪ-di -> Stress on 'ḥī'
  # 10. ʿalā (عَلَىٰ) -> ʿa-LĀ -> Stress on 'lā'
  
  stressed_syllables_list = [
      "نا",   # nā
      "بي",   # bay
      "ل",    # la
      "ضا",   # ḍan
      "في",   # fī
      "ذا",   # dhā
      "فص",   # faṣ
      "لآ",   # lā
      "حي",   # ḥī
      "لى"    # lā
  ]
  
  stressed_syllables_string = " ".join(stressed_syllables_list)

  # Print the final answer in the required format.
  print(author_name)
  print() # Leave exactly one line space
  print(stressed_syllables_string)

# Execute the function to generate the output.
get_solution()
<<<maimonides

نا بي ل ضا في ذا فص لآ حي لى>>>