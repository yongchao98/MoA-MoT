from itertools import permutations

# Books to arrange
books = ['A', 'B', 'C', 'E', 'F', 'G']

# Fixed position for Book D
fixed_positions = {2: 'D'}

# Function to check if a permutation satisfies all constraints
def is_valid(arrangement):
    # Place D in position 2
    arrangement = arrangement[:1] + ('D',) + arrangement[1:]
    
    # Constraint 1: Book A cannot be in position 6
    if arrangement[5] == 'A':
        return False
    
    # Constraint 2: Book E must be to the left of book C
    if arrangement.index('E') > arrangement.index('C'):
        return False
    
    # Constraint 3: Book F cannot be in position 4
    if arrangement[3] == 'F':
        return False
    
    # Constraint 4: Book C must be adjacent to book F
    c_index = arrangement.index('C')
    f_index = arrangement.index('F')
    if abs(c_index - f_index) != 1:
        return False
    
    # Constraint 5: Book A must be adjacent to book B
    a_index = arrangement.index('A')
    b_index = arrangement.index('B')
    if abs(a_index - b_index) != 1:
        return False
    
    # Constraint 6: Book E must be adjacent to book C
    e_index = arrangement.index('E')
    if abs(e_index - c_index) != 1:
        return False
    
    # Constraint 8: Book G must be to the left of book F
    if arrangement.index('G') > f_index:
        return False
    
    return True

# Generate all permutations of the books excluding D
for perm in permutations(books):
    if is_valid(perm):
        # Insert D in position 2
        result = list(perm[:1]) + ['D'] + list(perm[1:])
        print(result)
        break