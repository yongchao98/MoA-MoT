# Check if a modulo operation can relate the positions
# We will try a few modulo bases to see if any fit the pattern

# Known position
known_x, known_y = 1, 3

# Question mark position
question_x, question_y = 16, 10

# Try different modulo bases
for base in range(2, 20):
    if (known_x % base == question_x % base) and (known_y % base == question_y % base):
        print(f"Pattern found with modulo base: {base}")
        break
else:
    print("No simple modulo pattern found.")