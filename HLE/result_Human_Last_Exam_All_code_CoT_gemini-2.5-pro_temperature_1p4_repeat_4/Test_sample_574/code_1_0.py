import sys
import io

# Backup original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

# Step 1: Define the source words based on the puzzle's descriptions.
source_word_1 = "capes"
description_1 = "sleeveless garments that drape over the back and shoulders"

source_word_2 = "mentors"
description_2 = "experienced and trusted individuals who guide and advise others"

# Step 2: Form a word from the letters of the first source word.
# The letters of "capes" are c, a, p, e, s. We can form the word "space".
word_1 = "space"

# Step 3: Form a word from the letters of the second source word.
# The letters of "mentors" are m, e, n, t, o, r, s. We can form the word "storm".
word_2 = "storm"

# Step 4: Combine the two words to form the ship name starting with 'S'.
# The combination of "space" and "storm" creates "SpaceStorm".
ship_name = "SpaceStorm"

# Step 5: Print the solution, showing each component of the final "equation".
# This fulfills the prompt's request to show how the final answer is constructed.
print(f"The word from the letters of '{source_word_1}' ({description_1}) is: {word_1}")
print(f"The word from the letters of '{source_word_2}' ({description_2}) is: {word_2}")
print("\nCombining these gives the final ship name:")
print(f"{word_1} + {word_2} = {ship_name}")

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
output = captured_output.getvalue()
print(output)
