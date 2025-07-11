import sys
import io

# Ensure the output is printed in UTF-8 for Arabic characters.
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

# 1) The name of the author of the book.
author_name = "maimonides"

# 2) A list of the syllables with primary word stress for the first 10 words.
stressed_syllables = "بي ل أي في ها فص لا حي لى رأ"

# Print the final answer with the required formatting.
print(author_name)
print() # Exactly one line space between sections
print(stressed_syllables)