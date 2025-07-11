# The riddle refers to the Chinese character 福 (fú), which means good fortune.
# This character is composed of several smaller parts (radicals/components).
# This script represents the character's structure as an equation of the stroke counts
# of its main components, as per the user's request.

# Stroke counts for each component of the character 福
radical_shi_strokes = 4
component_one_strokes = 1
component_mouth_strokes = 3
component_field_strokes = 5

# Calculate the total number of strokes
total_strokes = radical_shi_strokes + component_one_strokes + component_mouth_strokes + component_field_strokes

# Print the final equation showing the sum of strokes
print(f"The character is 福 (fú). Its structure can be represented by the stroke counts of its parts:")
print(f"{radical_shi_strokes} + {component_one_strokes} + {component_mouth_strokes} + {component_field_strokes} = {total_strokes}")