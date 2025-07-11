# This script identifies two languages based on a set of linguistic clues.

# Clues for language 'a':
# 1. Is a recognized language still in use.
# 2. Orthography lacks the letters 'k' and 'w'.
# 3. Orthography includes the letter 'à'.
# Analysis: These clues point to Italian. Its standard 21-letter alphabet excludes 'k' and 'w',
# and it uses the grave accent on vowels, including 'à' (e.g., in "città").
language_a = "Italian"

# Clues for language 'b':
# 1. Is a recognized language still in use.
# 2. The letter combinations "ggj" and "skt" are very widely used.
# Analysis: These specific clusters are characteristic of Faroese.
# "ggj" is found in common words like "leggja" (to lay).
# "skt" is a common adjectival ending, as in "føroyskt" (Faroese).
language_b = "Faroese"

# Print the solution
print(f"Based on the analysis of the clues:")
print(f"Language a is: {language_a}")
print(f"Language b is: {language_b}")