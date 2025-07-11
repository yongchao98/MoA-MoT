# This script determines the two languages based on a set of orthographic and usage clues.

# Step 1: Identify Language 'a'.
# The clues are: no 'k' or 'w' in the orthography, but it contains 'à'.
# This points directly to Italian. The standard Italian alphabet has 21 letters,
# excluding 'k' and 'w'. The grave accent 'à' is commonly used (e.g., in 'città').
language_a = "Italian"

# Step 2: Identify Language 'b'.
# The clues are: the letter combinations "ggj" and "skt" are widely used.
# The "ggj" cluster is highly characteristic of Maltese, where it represents the
# common geminated (doubled) 'ġġ' sound (e.g., in 'oġġett').
# The "skt" cluster is also found in Maltese words like 'manuskritt' (manuscript).
language_b = "Maltese"

# Step 3: Print the results.
# The code below will display the identified languages.
print("Based on the clues provided:")
print(f"Language a is: {language_a}")
print(f"Language b is: {language_b}")