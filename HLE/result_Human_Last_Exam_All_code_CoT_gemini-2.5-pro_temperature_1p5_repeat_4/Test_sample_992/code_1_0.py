# The riddle describes the visual components of a Chinese character.
# "One horizontal stroke, another horizontal stroke..." suggests two horizontal strokes.
# "...one vertical on the left, one vertical on the right" describes two vertical strokes forming a frame.
# This structure perfectly matches the character '井' (jǐng), which means 'well'.
# The character '井' is composed of 4 strokes. We can model this by adding up the strokes mentioned in the riddle.

# Numbers representing the strokes mentioned in the riddle
h_stroke_1 = 1
h_stroke_2 = 1
v_stroke_1 = 1
v_stroke_2 = 1

# Calculate the total number of strokes
total = h_stroke_1 + h_stroke_2 + v_stroke_1 + v_stroke_2

print("The riddle describes the character '井' (jǐng).")
print("It has 2 horizontal and 2 vertical strokes, for a total of 4 strokes.")
print("We can represent this count with the following equation based on the riddle:")

# Print the full equation showing each number
print(f"Equation: {h_stroke_1} + {h_stroke_2} + {v_stroke_1} + {v_stroke_2} = {total}")

print("\nThe character is:")
print("井")
