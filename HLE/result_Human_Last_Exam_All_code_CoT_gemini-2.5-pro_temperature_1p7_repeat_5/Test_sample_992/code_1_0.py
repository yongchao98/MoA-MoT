# The riddle describes the Chinese character "王" (wáng), which means "king".
# The thinking process is to deconstruct the character into its constituent strokes.

# The first part of the riddle, "One horizontal stroke, another horizontal stroke, after another",
# clearly describes the three horizontal strokes (三).
horizontal_strokes = 3

# The rest of the riddle poetically describes the single vertical stroke
# that intersects the three horizontal ones.
vertical_strokes = 1

# The total number of strokes in the character "王" is the sum of its horizontal and vertical strokes.
total_strokes = horizontal_strokes + vertical_strokes

# Now, we print the explanation and the final equation as requested.
print("The character described is '王' (wáng), meaning 'king'.")
print("It is composed of the following strokes:")
print(f"Horizontal strokes: {horizontal_strokes}")
print(f"Vertical strokes: {vertical_strokes}")
print("\nThe equation representing the character's total stroke count is:")
print(f"{horizontal_strokes} + {vertical_strokes} = {total_strokes}")