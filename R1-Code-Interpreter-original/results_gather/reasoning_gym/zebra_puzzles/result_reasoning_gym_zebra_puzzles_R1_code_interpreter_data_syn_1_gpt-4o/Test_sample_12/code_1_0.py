from itertools import permutations

names = ['alice', 'arnold', 'bob', 'carol']
cigars = ['pall mall', 'dunhill', 'blue master', 'prince']
lunches = ['stir fry', 'grilled cheese', 'soup', 'pizza']
colors = ['blue', 'purple', 'brown', 'white']

# Generate all permutations of the lists
name_perms = permutations(names)
cigar_perms = permutations(cigars)
lunch_perms = permutations(lunches)
color_perms = permutations(colors)

# Check each permutation for the constraints
for name in name_perms:
    for cigar in cigar_perms:
        for lunch in lunch_perms:
            for color in color_perms:
                if (color[0] == 'white' and lunch[0] == 'stir fry' and cigar[0] == 'pall mall' and
                    color[2] == 'brown' and cigar[3] == 'prince' and
                    lunch[1] == 'pizza' and name[1] == 'arnold' and
                    name[2] == 'alice' and lunch[2] == 'soup' and color[2] == 'blue' and
                    cigar[2] == 'blue master' and cigar[3] == 'dunhill' and name[3] == 'bob'):
                    print(name[0])
                    break