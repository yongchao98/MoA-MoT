from itertools import permutations

# Fixed positions based on constraints
books = ['A', 'E', None, None, None, None, 'C']

# Possible books for positions 3, 4, 5, 6
remaining_books = ['D', 'G', 'B', 'F']

# Check all permutations of the remaining books
for perm in permutations(remaining_books):
    # Place the permutation in the available slots
    books[2], books[3], books[4], books[5] = perm
    
    # Check all constraints
    if books.index('D') < books.index('C') and \
       books.index('G') > books.index('D') and \
       books.index('B') > books.index('G'):
        print(books)
        break