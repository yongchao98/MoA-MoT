from itertools import permutations

# Define the possible attributes
names = ['arnold', 'bob', 'carol', 'alice']
nationalities = ['brit', 'chinese', 'dane', 'german']
colors = ['purple', 'white', 'blue', 'brown']
drinks = ['milk', 'water', 'boba tea', 'coffee']

# Define the constraints based on the clues
def is_valid(arrangement):
    # Unpack the arrangement
    (n1, n2, n3, n4), (nat1, nat2, nat3, nat4), (c1, c2, c3, c4), (d1, d2, d3, d4) = arrangement
    
    # Apply the constraints
    if n1 != 'arnold' and n2 != 'arnold' and n3 != 'arnold':
        return False
    if n1 == 'arnold' and c1 != 'purple':
        return False
    if n2 == 'arnold' and c2 != 'purple':
        return False
    if n3 == 'arnold' and c3 != 'purple':
        return False
    if n4 == 'arnold':
        return False
    if n1 == 'arnold' and d1 != 'coffee':
        return False
    if n2 == 'arnold' and d2 != 'coffee':
        return False
    if n3 == 'arnold' and d3 != 'coffee':
        return False
    if n4 == 'arnold':
        return False
    if n1 == 'bob' and d1 != 'boba tea':
        return False
    if n2 == 'bob' and d2 != 'boba tea':
        return False
    if n3 == 'bob' and d3 != 'boba tea':
        return False
    if n4 == 'bob':
        return False
    if n1 == 'carol' and c1 != 'white':
        return False
    if n2 == 'carol' and c2 != 'white':
        return False
    if n3 == 'carol' and c3 != 'white':
        return False
    if n4 == 'carol' and c4 != 'white':
        return False
    if n1 == 'carol' and nat1 != 'chinese':
        return False
    if n2 == 'carol' and nat2 != 'chinese':
        return False
    if n3 == 'carol' and nat3 != 'chinese':
        return False
    if n4 == 'carol' and nat4 != 'chinese':
        return False
    if nat1 == 'german' and nat2 != 'dane':
        return False
    if nat2 == 'german' and nat3 != 'dane':
        return False
    if nat3 == 'german' and nat4 != 'dane':
        return False
    if nat4 == 'german':
        return False
    if nat1 == 'brit' and c1 != 'brown':
        return False
    if nat2 == 'brit' and c2 != 'brown':
        return False
    if nat3 == 'brit' and c3 != 'brown':
        return False
    if nat4 == 'brit' and c4 != 'brown':
        return False
    if c1 == 'purple' and d2 != 'milk':
        return False
    if c2 == 'purple' and d3 != 'milk':
        return False
    if c3 == 'purple' and d4 != 'milk':
        return False
    if c4 == 'purple':
        return False
    if n1 == 'bob' and n2 != 'alice':
        return False
    if n2 == 'bob' and n3 != 'alice':
        return False
    if n3 == 'bob' and n4 != 'alice':
        return False
    if n4 == 'bob':
        return False
    
    return True

# Generate all possible permutations of the attributes
for arrangement in permutations(names):
    for nationality in permutations(nationalities):
        for color in permutations(colors):
            for drink in permutations(drinks):
                if is_valid((arrangement, nationality, color, drink)):
                    print(arrangement)
                    break