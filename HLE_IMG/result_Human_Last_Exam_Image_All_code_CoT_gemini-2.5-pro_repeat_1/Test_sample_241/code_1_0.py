import sys
import io

# Set stdout to utf-8 to handle Arabic characters
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

def get_solution():
  """
  This function contains the solution to the problem.
  Part 1: Identifies the author of the text.
  Part 2: Lists the stressed syllables of the first 10 words.
  """
  
  # Part 1: The author's name.
  # The text is from Chapter 75 of Maimonides' "Guide for the Perplexed".
  # The author is Moses Maimonides.
  author = "maimonides"
  
  # Part 2: Stressed syllables.
  # The first 10 words in Arabic are:
  # أنا أبين لك أيضا في هذا الفصل ودلائل التوحيد على
  # (anā ubayyinu laka ayḍan fī hādhā al-faṣli wa-dalāʾili al-tawḥīdi ʿalā)
  #
  # The stressed syllables (unvocalized) are:
  # 1. anā -> nā -> نا
  # 2. ubayyinu -> bay -> بي
  # 3. laka -> la -> ل
  # 4. ayḍan -> ay -> اي
  # 5. fī -> fī -> في
  # 6. hādhā -> hā -> ها
  # 7. al-faṣli -> faṣ -> فص
  # 8. dalāʾili -> lā -> لا
  # 9. al-tawḥīdi -> ḥī -> حي
  # 10. ʿalā -> lā -> لى
  stressed_syllables = "نا بي ل اي في ها فص لا حي لى"
  
  # Printing the final answer as requested
  print(author)
  print() # Exactly one line space
  print(stressed_syllables)

get_solution()
