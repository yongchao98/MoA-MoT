# The riddle describes the creation of a Chinese character.

# "One horizontal stroke, another horizontal stroke..."
horizontal_strokes = 2

# "...one vertical on the left, one vertical on the right"
vertical_strokes = 2

# The riddle asks for the character formed by these strokes.
# To fulfill the requirement of showing an equation, we can sum the strokes.
total_strokes = horizontal_strokes + vertical_strokes

# Print each number in the final equation.
print(f"{horizontal_strokes} + {vertical_strokes} = {total_strokes}")

# The character with two horizontal and two vertical strokes is "井".
print("The character is: 井")