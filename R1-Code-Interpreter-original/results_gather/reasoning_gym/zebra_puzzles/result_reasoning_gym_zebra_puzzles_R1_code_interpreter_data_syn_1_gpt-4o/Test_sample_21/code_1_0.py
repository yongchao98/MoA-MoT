from itertools import permutations

names = ['carol', 'bob', 'alice', 'arnold']
colors = ['brown', 'white', 'blue', 'purple']
nationalities = ['dane', 'brit', 'german', 'chinese']
flowers = ['daffodils', 'carnations', 'lilies', 'iris']

# Constraints
def is_valid(arrangement):
    # Unpack the arrangement
    (name1, name2, name3, name4) = arrangement['names']
    (color1, color2, color3, color4) = arrangement['colors']
    (nationality1, nationality2, nationality3, nationality4) = arrangement['nationalities']
    (flower1, flower2, flower3, flower4) = arrangement['flowers']
    
    # Apply the constraints
    if nationality3 != 'german' or color3 != 'purple':
        return False
    if name2 != 'carol' or flower2 != 'carnations':
        return False
    if not ((name2 == 'carol' and nationality1 == 'brit') or (name3 == 'carol' and nationality4 == 'brit')):
        return False
    if not ((name1 == 'alice' and nationality2 == 'chinese') or (name2 == 'alice' and nationality3 == 'chinese') or (name3 == 'alice' and nationality4 == 'chinese')):
        return False
    if name3 != 'arnold' or color3 != 'blue':
        return False
    if flower4 != 'daffodils':
        return False
    if nationality1 != 'dane' or name1 != 'bob':
        return False
    if not ((color1 == 'brown' and flower2 == 'lilies') or (color2 == 'brown' and flower3 == 'lilies') or (color3 == 'brown' and flower4 == 'lilies')):
        return False
    return True

# Try all permutations
for names_perm in permutations(names):
    for colors_perm in permutations(colors):
        for nationalities_perm in permutations(nationalities):
            for flowers_perm in permutations(flowers):
                arrangement = {
                    'names': names_perm,
                    'colors': colors_perm,
                    'nationalities': nationalities_perm,
                    'flowers': flowers_perm
                }
                if is_valid(arrangement):
                    print(arrangement)
                    break