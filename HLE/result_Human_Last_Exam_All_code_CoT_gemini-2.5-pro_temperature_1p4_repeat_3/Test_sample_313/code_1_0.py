# This script identifies two languages, 'a' and 'b', based on a set of linguistic clues.

# --- Language 'a' Identification ---
# The clues for language 'a' are:
# 1. It is a living language, recognized officially or as a minority language.
# 2. Its orthography does not contain the letters 'k' or 'w'.
# 3. Its orthography contains the letter 'à'.
#
# Analysis: The absence of 'k' and 'w' from the native alphabet is a strong clue.
# Italian fits perfectly: its standard alphabet consists of 21 letters and does not
# include 'j', 'k', 'w', 'x', or 'y'. The grave accent 'à' is commonly used
# in Italian to mark stress on the final vowel (e.g., 'città').
# Therefore, language 'a' is identified as Italian.
language_a = "Italian"


# --- Language 'b' Identification ---
# The clues for language 'b' are:
# 1. It is a living language, recognized officially or as a minority language.
# 2. The letter combination "ggj" is very widely used.
# 3. The letter combination "skt" is very widely used.
#
# Analysis: The letter clusters "ggj" and "skt" are highly distinctive.
# These are characteristic of the Faroese language.
# - "ggj" appears in many common verbs like 'byggja' (to build) and 'leggja' (to lay).
# - "-skt" is the standard neuter ending for adjectives, including names of languages
#   like 'føroyskt' (Faroese) and 'enskt' (English).
# Therefore, language 'b' is identified as Faroese.
language_b = "Faroese"

# --- Final Answer ---
print("Based on the clues:")
print("Language a is: " + language_a)
print("Language b is: " + language_b)