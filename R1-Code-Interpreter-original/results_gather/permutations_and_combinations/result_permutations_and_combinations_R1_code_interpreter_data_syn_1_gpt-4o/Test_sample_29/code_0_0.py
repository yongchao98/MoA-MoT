# Initial arrangement based on constraints
arrangement = ["C", "D", "A", "E", None, None]

# Remaining books
remaining_books = ["B", "F"]

# Fill the remaining positions
arrangement[4] = remaining_books[0]
arrangement[5] = remaining_books[1]

print(arrangement)