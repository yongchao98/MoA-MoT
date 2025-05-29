from itertools import permutations

# Define the possible attributes
names = ['alice', 'arnold', 'bob', 'carol']
cigars = ['pall mall', 'dunhill', 'blue master', 'prince']
lunches = ['stir fry', 'grilled cheese', 'soup', 'pizza']
colors = ['blue', 'purple', 'brown', 'white']

# Define the constraints
def is_valid(arrangement):
    # Unpack the arrangement
    (name1, name2, name3, name4) = arrangement['names']
    (cigar1, cigar2, cigar3, cigar4) = arrangement['cigars']
    (lunch1, lunch2, lunch3, lunch4) = arrangement['lunches']
    (color1, color2, color3, color4) = arrangement['colors']
    
    # Apply the constraints
    if name1 == 'bob' and cigar1 != 'dunhill':
        return False
    if name2 == 'bob' and cigar2 != 'dunhill':
        return False
    if name3 == 'bob' and cigar3 != 'dunhill':
        return False
    if name4 == 'bob' and cigar4 != 'dunhill':
        return False
    
    if name1 == 'alice' and lunch1 != 'soup':
        return False
    if name2 == 'alice' and lunch2 != 'soup':
        return False
    if name3 == 'alice' and lunch3 != 'soup':
        return False
    if name4 == 'alice' and lunch4 != 'soup':
        return False
    
    if name1 == 'alice' and color1 != 'blue':
        return False
    if name2 == 'alice' and color2 != 'blue':
        return False
    if name3 == 'alice' and color3 != 'blue':
        return False
    if name4 == 'alice' and color4 != 'blue':
        return False
    
    if cigar1 == 'pall mall' and color1 != 'white':
        return False
    if cigar2 == 'pall mall' and color2 != 'white':
        return False
    if cigar3 == 'pall mall' and color3 != 'white':
        return False
    if cigar4 == 'pall mall' and color4 != 'white':
        return False
    
    if cigar1 == 'pall mall' and lunch1 != 'stir fry':
        return False
    if cigar2 == 'pall mall' and lunch2 != 'stir fry':
        return False
    if cigar3 == 'pall mall' and lunch3 != 'stir fry':
        return False
    if cigar4 == 'pall mall' and lunch4 != 'stir fry':
        return False
    
    if color1 == 'brown' and cigar2 != 'prince':
        return False
    if color2 == 'brown' and cigar3 != 'prince':
        return False
    if color3 == 'brown' and cigar4 != 'prince':
        return False
    
    if lunch1 == 'pizza' and name2 != 'arnold' and name1 != 'arnold':
        return False
    if lunch2 == 'pizza' and name3 != 'arnold' and name2 != 'arnold':
        return False
    if lunch3 == 'pizza' and name4 != 'arnold' and name3 != 'arnold':
        return False
    
    if cigar1 == 'blue master' and cigar2 != 'dunhill':
        return False
    if cigar2 == 'blue master' and cigar3 != 'dunhill':
        return False
    if cigar3 == 'blue master' and cigar4 != 'dunhill':
        return False
    
    return True

# Generate all possible arrangements
for name_perm in permutations(names):
    for cigar_perm in permutations(cigars):
        for lunch_perm in permutations(lunches):
            for color_perm in permutations(colors):
                arrangement = {
                    'names': name_perm,
                    'cigars': cigar_perm,
                    'lunches': lunch_perm,
                    'colors': color_perm
                }
                if is_valid(arrangement):
                    print(arrangement)
                    break