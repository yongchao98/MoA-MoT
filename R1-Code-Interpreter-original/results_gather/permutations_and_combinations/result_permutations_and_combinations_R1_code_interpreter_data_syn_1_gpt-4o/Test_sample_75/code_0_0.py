from itertools import permutations

# Positions are 1 to 7
positions = [1, 2, 3, 4, 5, 6, 7]

# Fixed positions based on constraints
fixed_positions = {
    'D': 5,
    'B': 4,
    'C': 6,
    'E': 7
}

# Remaining books to place
remaining_books = ['A', 'G', 'F']

# Generate permutations for the remaining books
for perm in permutations(remaining_books):
    # Create a full arrangement
    arrangement = [''] * 7
    arrangement[fixed_positions['D'] - 1] = 'D'
    arrangement[fixed_positions['B'] - 1] = 'B'
    arrangement[fixed_positions['C'] - 1] = 'C'
    arrangement[fixed_positions['E'] - 1] = 'E'
    
    # Place the remaining books
    j = 0
    for i in range(7):
        if arrangement[i] == '':
            arrangement[i] = perm[j]
            j += 1
    
    # Check constraints
    if arrangement.index('A') < arrangement.index('D') and arrangement.index('G') < arrangement.index('F'):
        print(arrangement)
        break