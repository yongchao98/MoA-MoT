from itertools import permutations

# Books to arrange
books = ['B', 'C', 'D', 'E', 'F', 'G']

# Constraint: A is in position 7
# We need to find a valid permutation of the remaining books in positions 1 to 6
for perm in permutations(books):
    # Check all constraints
    G_index = perm.index('G')
    C_index = perm.index('C')
    E_index = perm.index('E')
    F_index = perm.index('F')
    D_index = perm.index('D')
    B_index = perm.index('B')
    
    if (G_index != 2 and  # G cannot be in position 3 (index 2)
        F_index > E_index and  # F must be to the right of E
        D_index > F_index and  # D must be to the right of F
        C_index < G_index and  # C must be to the left of G
        E_index > B_index):  # E must be to the right of B
        # If all constraints are satisfied, we have found a valid arrangement
        result = list(perm) + ['A']  # Add A at the end
        print(result)
        break