# Step 1: Define the source words based on the riddle's clues.
garments_word = "CAPES"
individuals_word = "THEOLOGIANS"

# Step 2: Form the two parts of the ship's name from the letters of the source words.
# "SCAPE" is formed from the letters s, c, a, p, e in "CAPES".
part1 = "SCAPE"

# "GOAT" is formed from the letters g, o, a, t in "THEOLOGIANS".
part2 = "GOAT"

# Step 3: Combine the two parts to get the final ship name.
ship_name = part1 + part2

# Step 4: Print the solution, showing how the parts combine to form the final name.
# The riddle asks to show each part of the final "equation".
print(f"The first word, derived from '{garments_word}', is '{part1}'.")
print(f"The second word, derived from '{individuals_word}', is '{part2}'.")
print("Combining them gives the Culture ship name:")
print(f"{part1} + {part2} = {ship_name}")