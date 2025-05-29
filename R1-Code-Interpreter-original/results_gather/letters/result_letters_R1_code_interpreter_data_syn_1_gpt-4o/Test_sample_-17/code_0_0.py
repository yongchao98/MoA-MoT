# Given string
word = 'bmmgtjgzorlbgkwvroerqhkg'

# Initialize count and positions list
count = 0
positions = []

# Iterate through the string
for index, char in enumerate(word):
    if char == 'g':
        count += 1
        positions.append(index + 1)  # Convert 0-based index to 1-based position

# Print the result
print(f"Count: {count}, Positions: {positions}")