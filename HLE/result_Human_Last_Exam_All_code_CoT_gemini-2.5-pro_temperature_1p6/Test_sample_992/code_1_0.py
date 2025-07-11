# The riddle describes the components of a Chinese character based on its strokes.

# "One horizontal stroke, another horizontal stroke, after another"
# This describes two horizontal strokes.
num_horizontal_strokes = 2

# "one vertical on the left, one vertical on the right"
# This describes two vertical strokes.
num_vertical_strokes = 2

# The character '井' (jǐng) is composed of these four strokes.
# We can represent this composition as a simple equation of the stroke counts.
total_strokes = num_horizontal_strokes + num_vertical_strokes

# Print the analysis and the equation of the stroke count.
print("The riddle describes a character composed of:")
print(f"- {num_horizontal_strokes} horizontal strokes")
print(f"- {num_vertical_strokes} vertical strokes")

print("\nThe stroke count can be represented by the equation:")
print(f"'{num_horizontal_strokes} horizontal strokes' + '{num_vertical_strokes} vertical strokes' = '{total_strokes} total strokes'")

# Print each number in the final equation as requested.
print(f"The numbers in this equation are: {num_horizontal_strokes}, {num_vertical_strokes}, and {total_strokes}.")

print("\nThe Chinese character that matches this description is:")
print("井")