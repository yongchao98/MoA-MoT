import sys
import io

# Set stdout to utf-8 to ensure Arabic characters are printed correctly.
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

def get_solution():
  """
  This function provides the solution based on the analysis of the Judaeo-Arabic text.
  """

  # Part 1: The name of the author of the book.
  # The text is identified as Chapter 75 of Maimonides' "Guide for the Perplexed".
  author_name = "maimonides"

  # Part 2: A list of the syllables with primary word stress for the first 10 words.
  # The first 10 words are: أنا أبيّن لك أيضاً في هذا الفصل دلائل التوحيد على
  # anā, ubayyinu, laka, ayḍan, fī, hādhā, al-faṣl(i), dalāʾil, al-tawḥīd(i), ʿalā
  # The stressed syllables are: nā, bay, la, ay, fī, hā, faṣ, lā, ḥī, lā
  stressed_syllables = "نا بي ل أي في ها فص لا حي لى"

  # Printing the final answer in the required format.
  print(author_name)
  print()  # For the single empty line
  print(stressed_syllables)

get_solution()
# The final answer is then wrapped in <<<>>>
final_answer_string = "maimonides\n\nنا بي ل أي في ها فص لا حي لى"
# print(f"<<<{final_answer_string}>>>")