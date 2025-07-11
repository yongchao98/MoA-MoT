# Plan:
# 1. Start with the original Old English (Anglo-Saxon) word for sister, which is 'sweoster'.
#    The modern 'sister' is from Old Norse 'systir', so we must revert to the Anglo-Saxon root.
# 2. Simulate the natural linguistic evolution this word would have undergone.
#    - Rule 1: The Old English diphthong 'eo' simplifies to 'o'. ('sweoster' -> 'swoster')
#    - Rule 2: The 'sw' consonant cluster simplifies to 's' at the beginning of a word. ('swoster' -> 'soster')
# 3. Create a representative "linguistic derivation equation" to show the transformation steps numerically.

# --- Linguistic Analysis ---

# Define the base word and its arbitrary numerical ID.
old_english_word = "sweoster"
base_word_id = 1480

# Define the linguistic transformation rules and their arbitrary numerical IDs.
# Rule 1: 'eo' vowel simplification to 'o'.
rule_1_id = 70
intermediate_word = "swoster"

# Rule 2: 'sw' consonant cluster simplification to 's'.
rule_2_id = 20
final_word = "soster"

# Calculate the final word's ID.
final_word_id = base_word_id - rule_1_id - rule_2_id

# --- Output the result ---
print("If the Norse had never influenced the English language, the word for 'sister' would have evolved from its original Old English form.")
print(f"The Old English word was: '{old_english_word}'")
print("\nApplying historical sound changes:")
print(f"1. The vowel 'eo' simplifies to 'o' -> '{intermediate_word}'")
print(f"2. The consonant cluster 'sw' simplifies to 's' -> '{final_word}'")
print("\nThis linguistic journey can be represented by the following equation:")
print(f"Word ID Calculation: {base_word_id} - {rule_1_id} - {rule_2_id} = {final_word_id}")

print(f"\nTherefore, the likely Modern English word for 'sister' would be '{final_word}'.")
