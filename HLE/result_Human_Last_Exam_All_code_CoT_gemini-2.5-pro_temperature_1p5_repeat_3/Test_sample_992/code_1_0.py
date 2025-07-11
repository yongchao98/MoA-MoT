# The riddle describes how to form a Chinese character by its strokes.

# Step 1: Decode the first clue.
# "One horizontal stroke, another horizontal stroke, after another" describes
# the character for the number 3, which is composed of three horizontal strokes.
first_component = "三"

# Step 2: Decode the second, more complex clue.
# The rest of the riddle describes adding a single vertical stroke.
# This stroke intersects the three horizontal ones.
second_component = "丨"

# Step 3: Combine the components.
# Adding a vertical stroke through the middle of the character for three
# creates the character for "king".
final_character = "王"

# Step 4: Print the thinking process and the final equation.
print(f"'{first_component}' (three) is described by 'One horizontal stroke, another horizontal stroke, after another.'")
print(f"'{second_component}' (a vertical stroke) is poetically described by the rest of the riddle.")
print("When you combine them, you get the answer:")
print(f"{first_component} + {second_component} = {final_character}")
