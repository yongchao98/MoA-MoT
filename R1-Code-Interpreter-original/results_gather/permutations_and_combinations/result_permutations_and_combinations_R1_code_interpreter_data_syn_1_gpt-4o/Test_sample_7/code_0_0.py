import itertools

books = ['A', 'B', 'C', 'D', 'E', 'F']
constraints = [
    lambda p: p.index('D') > p.index('C'),  # D is to the right of C
    lambda p: p.index('F') < p.index('A'),  # F is to the left of A
    lambda p: p.index('F') < p.index('B'),  # F is to the left of B
    lambda p: p.index('C') != 3,            # C is not in position 4
    lambda p: p.index('A') != 2,            # A is not in position 3
    lambda p: p.index('E') != 2,            # E is not in position 3
    lambda p: p.index('C') != 1,            # C is not in position 2
    lambda p: p.index('F') < p.index('D')   # F is to the left of D
]

def is_valid_permutation(p):
    return all(constraint(p) for constraint in constraints)

for permutation in itertools.permutations(books):
    if is_valid_permutation(permutation):
        print(list(permutation))
        break