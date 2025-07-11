# Step 1: Define the source word for the first part of the riddle.
source_word_1_clue = "sleeveless garments that drape over the back and shoulders"
source_word_1 = "CAPES"

# Step 2: Form the first word from the letters of the source word.
derived_word_1 = "SCAPE"

# Step 3: Define the source phrase for the second part of the riddle.
source_word_2_clue = "experienced and trusted individuals who guide and advise others"
source_word_2 = "CARE GIVERS"

# Step 4: Form the second word from the letters of the source phrase.
derived_word_2 = "GRACE"

# Step 5: Combine the two derived words to get the final ship name.
ship_name = derived_word_1 + derived_word_2

# Print the step-by-step solution
print(f"The word for '{source_word_1_clue}' is '{source_word_1}'.")
print(f"From the letters of '{source_word_1}', we can form the word '{derived_word_1}'.")
print("-" * 20)
print(f"The words for '{source_word_2_clue}' are '{source_word_2}'.")
print(f"From the letters of '{source_word_2}', we can form the word '{derived_word_2}'.")
print("-" * 20)
print("Combining these two words gives the Culture series ship name.")
# The prompt requests to output each part of the final equation.
print(f"Final Equation: {derived_word_1} + {derived_word_2} = {ship_name}")