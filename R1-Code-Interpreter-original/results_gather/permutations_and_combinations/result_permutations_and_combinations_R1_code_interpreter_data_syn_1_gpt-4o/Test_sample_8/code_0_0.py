from itertools import permutations

books = ['A', 'B', 'C', 'D', 'E', 'F']

# F is fixed at position 6
for perm in permutations(books[:-1]):  # permute A, B, C, D, E
    if perm[2] == 'C':  # C cannot be in position 3
        continue
    if perm[2] == 'A':  # A cannot be in position 3
        continue
    if 'A' in perm and 'D' in perm:
        a_index = perm.index('A')
        d_index = perm.index('D')
        if abs(a_index - d_index) != 1:  # A and D must be adjacent
            continue
    if 'C' in perm and 'F' in perm:
        c_index = perm.index('C')
        if c_index > 4:  # C must be before F
            continue
    if 'A' in perm and 'F' in perm:
        a_index = perm.index('A')
        if a_index > 4:  # A must be before F
            continue
    if 'B' in perm and 'E' in perm:
        b_index = perm.index('B')
        e_index = perm.index('E')
        if b_index > e_index:  # B must be before E
            continue
    # If all conditions are satisfied, we have a valid arrangement
    result = list(perm) + ['F']
    print(result)
    break