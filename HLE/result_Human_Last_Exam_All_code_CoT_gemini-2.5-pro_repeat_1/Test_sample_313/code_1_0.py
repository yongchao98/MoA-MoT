# This script identifies the two languages based on the provided clues.

# Step 1: Identify Language 'a'.
# Clues: No 'k' or 'w', but has 'à'. Still in use.
# Analysis: The standard Italian alphabet for native words omits 'k' and 'w'.
# It frequently uses 'à' to mark stress on the final vowel (e.g., "città").
# Therefore, Language 'a' is identified as Italian.
language_a = "Italian"

# Step 2: Identify Language 'b'.
# Clues: Widely uses "ggj" and "skt". Still in use.
# Analysis: The letter combination "ggj" is a distinctive feature of Faroese
# (e.g., "leggja"). The cluster "skt" is also common, particularly in
# nationality adjectives (e.g., "føroyskt").
# Therefore, Language 'b' is identified as Faroese.
language_b = "Faroese"

# Step 3: Print the final answer.
print("Based on the linguistic clues provided:")
print(f"Language a is: {language_a}")
print(f"Language b is: {language_b}")