from itertools import permutations

# Possible names, colors, children, and smoothies
names = ['carol', 'arnold', 'alice', 'bob']
colors = ['white', 'purple', 'blue', 'brown']
children = ['timothy', 'bella', 'alice', 'billy']
smoothies = ['darkness', 'dragonfruit', 'desert', 'butterscotch']

# Constraints based on the clues
def is_valid(arrangement):
    # Unpack the arrangement
    (name1, name2, name3, name4) = arrangement['names']
    (color1, color2, color3, color4) = arrangement['colors']
    (child1, child2, child3, child4) = arrangement['children']
    (smoothie1, smoothie2, smoothie3, smoothie4) = arrangement['smoothies']
    
    # Apply the clues
    if child3 != 'alice' or name4 != 'arnold':  # Clue 1
        return False
    if smoothie1 != 'darkness' or color1 != 'white':  # Clue 2
        return False
    if child1 != 'timothy' or color2 != 'brown':  # Clue 3
        return False
    if name1 != 'bob' or color1 != 'white':  # Clue 4
        return False
    if child2 != 'bella':  # Clue 5
        return False
    if smoothie1 != 'desert' or child2 != 'bella':  # Clue 6
        return False
    if color1 != 'blue' or color2 != 'purple':  # Clue 7
        return False
    if name3 != 'carol' or smoothie3 != 'dragonfruit':  # Clue 8
        return False
    
    return True

# Generate all permutations and check for a valid arrangement
for name_perm in permutations(names):
    for color_perm in permutations(colors):
        for child_perm in permutations(children):
            for smoothie_perm in permutations(smoothies):
                arrangement = {
                    'names': name_perm,
                    'colors': color_perm,
                    'children': child_perm,
                    'smoothies': smoothie_perm
                }
                if is_valid(arrangement):
                    print(arrangement)
                    break