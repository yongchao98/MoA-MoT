from itertools import permutations

names = ['bob', 'alice', 'carol', 'arnold']
colors = ['blue', 'brown', 'white', 'purple']
smoothies = ['desert', 'darkness', 'butterscotch', 'dragonfruit']
children = ['billy', 'alice', 'timothy', 'bella']

# Constraints
def is_valid(arrangement):
    # Unpack the arrangement
    (name1, name2, name3, name4) = arrangement['names']
    (color1, color2, color3, color4) = arrangement['colors']
    (smoothie1, smoothie2, smoothie3, smoothie4) = arrangement['smoothies']
    (child1, child2, child3, child4) = arrangement['children']
    
    # Apply the constraints
    if name4 != 'carol':
        return False
    if color1 == 'white' and child1 != 'alice':
        return False
    if color2 == 'white' and child2 != 'alice':
        return False
    if color3 == 'white' and child3 != 'alice':
        return False
    if color4 == 'white' and child4 != 'alice':
        return False
    if name1 == 'bob' and color1 != 'purple':
        return False
    if name2 == 'bob' and color2 != 'purple':
        return False
    if name3 == 'bob' and color3 != 'purple':
        return False
    if name4 == 'bob' and color4 != 'purple':
        return False
    if name1 == 'arnold' and (color1 != 'blue' or smoothie1 != 'dragonfruit'):
        return False
    if name2 == 'arnold' and (color2 != 'blue' or smoothie2 != 'dragonfruit'):
        return False
    if name3 == 'arnold' and (color3 != 'blue' or smoothie3 != 'dragonfruit'):
        return False
    if name4 == 'arnold' and (color4 != 'blue' or smoothie4 != 'dragonfruit'):
        return False
    if smoothie1 == 'butterscotch' and child1 != 'timothy':
        return False
    if smoothie2 == 'butterscotch' and child2 != 'timothy':
        return False
    if smoothie3 == 'butterscotch' and child3 != 'timothy':
        return False
    if smoothie4 == 'butterscotch' and child4 != 'timothy':
        return False
    if child2 != 'billy':
        return False
    if color1 == 'brown' and smoothie2 != 'butterscotch':
        return False
    if color2 == 'brown' and smoothie3 != 'butterscotch':
        return False
    if color3 == 'brown' and smoothie4 != 'butterscotch':
        return False
    if color4 == 'brown':
        return False
    if color1 == 'white' and smoothie1 != 'desert':
        return False
    if color2 == 'white' and smoothie2 != 'desert':
        return False
    if color3 == 'white' and smoothie3 != 'desert':
        return False
    if color4 == 'white' and smoothie4 != 'desert':
        return False
    return True

# Generate all permutations and check constraints
for name_perm in permutations(names):
    for color_perm in permutations(colors):
        for smoothie_perm in permutations(smoothies):
            for child_perm in permutations(children):
                arrangement = {
                    'names': name_perm,
                    'colors': color_perm,
                    'smoothies': smoothie_perm,
                    'children': child_perm
                }
                if is_valid(arrangement):
                    print(arrangement)
                    break