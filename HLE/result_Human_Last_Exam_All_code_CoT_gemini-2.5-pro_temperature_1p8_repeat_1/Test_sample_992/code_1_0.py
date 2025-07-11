# This program solves a Chinese character riddle by breaking down the character's components.
# The riddle describes a fence, which points to the character '栅' (zhà).

# The character '栅' is composed of a left part and a right part.
# Left Part: 木 (mù), meaning 'wood'.
# Right Part: 册 (cè), meaning 'booklet', visually resembling fence slats.

# Let's define the number of strokes for each part.
strokes_in_wood_radical = 4
strokes_in_booklet_radical = 5

# The total number of strokes is the sum of the strokes in its components.
total_strokes = strokes_in_wood_radical + strokes_in_booklet_radical

print("The riddle describes the character '栅' (zhà), which means 'fence'.")
print("It's made of a 'wood' radical (木) on the left and a 'booklet' radical (册) on the right.")
print("The total number of strokes is the sum of the strokes of its two parts.")
print("Final Equation:")

# Print each number and symbol of the final equation as requested.
print(strokes_in_wood_radical, '+', strokes_in_booklet_radical, '=', total_strokes)

print("\nThe character is: 栅")