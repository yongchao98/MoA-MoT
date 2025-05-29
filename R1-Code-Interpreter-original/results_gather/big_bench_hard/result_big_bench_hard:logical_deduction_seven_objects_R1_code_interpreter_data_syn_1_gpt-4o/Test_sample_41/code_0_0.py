# Initialize the positions of the books
books = [''] * 7

# Place the books based on the given conditions
books[1] = 'Red'    # Second from the left
books[2] = 'Yellow' # Third from the left
books[3] = 'White'  # Fourth from the left
books[5] = 'Black'  # Second from the right
books[4] = 'Gray'   # Third from the right
books[6] = 'Green'  # Must be to the right of Red

# The only position left for the Purple book is the first position
books[0] = 'Purple'

# Print the arrangement of books
print(books)