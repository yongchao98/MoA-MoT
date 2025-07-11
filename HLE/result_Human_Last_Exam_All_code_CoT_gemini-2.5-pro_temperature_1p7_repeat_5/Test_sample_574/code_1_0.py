# 1. Define the source and derived words based on the riddle's solution.
source_garment = "SLEEPWEAR"
word1 = "SLEEPER"

source_guide = "RECEIVERS"
word2 = "SERVICE"

# 2. Combine the derived words to form the final ship name.
# The ship name is a concatenation of the two derived words.
ship_name = f"{word1} {word2}"

# 3. Print the solution's logic step-by-step.
print(f"The word for 'sleeveless garments' is '{source_garment}'.")
print(f"From its letters, we can make the word: '{word1}'")
print("-" * 30)
print(f"The word for 'trusted individuals who guide and advise' is '{source_guide}'.")
print(f"From its letters, we can make the word: '{word2}'")
print("-" * 30)
print("These two words are combined to form the Culture series ship name.")

# 4. Print the final equation, showing each word component.
print("\nThe final equation is:")
print(f"{word1} + {word2} = {ship_name}")