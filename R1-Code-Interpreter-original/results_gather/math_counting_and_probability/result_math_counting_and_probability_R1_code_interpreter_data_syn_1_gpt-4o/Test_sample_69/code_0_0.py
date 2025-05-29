# Total number of books
total_books = 6

# Calculate the number of valid distributions
# We exclude the cases where all books are in one location
valid_distributions = total_books - 1  # Excluding all in library and all checked out

print(valid_distributions)