import math

# Calculate the square root of 44/95
value = math.sqrt(44 / 95)

# Prepare the results for both positive and negative roots
positive_y = value
negative_y = -value

# Round the results to 4 decimal places
positive_y_rounded = round(positive_y, 4)
negative_y_rounded = round(negative_y, 4)

# Print the results
print(negative_y_rounded, positive_y_rounded)