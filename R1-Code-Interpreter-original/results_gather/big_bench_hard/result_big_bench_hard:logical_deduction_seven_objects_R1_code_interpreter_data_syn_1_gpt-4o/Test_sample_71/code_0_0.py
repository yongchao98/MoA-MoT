from itertools import permutations

# Fixed positions
positions = ['White', None, None, 'Gray', None, 'Orange', 'Black']

# Books to place
books = ['Red', 'Blue', 'Brown']

# Find the correct permutation
for perm in permutations(books):
    # Place the books in the available positions
    temp_positions = positions[:]
    temp_positions[1] = perm[0]
    temp_positions[2] = perm[1]
    temp_positions[4] = perm[2]
    
    # Check conditions
    if temp_positions.index('Red') > temp_positions.index('Blue') and temp_positions.index('Gray') > temp_positions.index('Red'):
        print(temp_positions)
        break