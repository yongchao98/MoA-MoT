from itertools import permutations

# Define the possible options
names = ['arnold', 'alice', 'bob', 'carol']
children = ['alice', 'bella', 'timothy', 'billy']
animals = ['bird', 'cat', 'fish', 'dog']
smoothies = ['desert', 'butterscotch', 'dragonfruit', 'darkness']

# Define the constraints
def is_valid(arrangement):
    # Unpack the arrangement
    (name1, name2, name3, name4) = arrangement['names']
    (child1, child2, child3, child4) = arrangement['children']
    (animal1, animal2, animal3, animal4) = arrangement['animals']
    (smoothie1, smoothie2, smoothie3, smoothie4) = arrangement['smoothies']
    
    # Apply the constraints
    if name3 != 'arnold': return False
    if 'bob' not in arrangement['names']: return False
    if smoothie1 == 'darkness' or smoothie2 == 'darkness' or smoothie3 == 'darkness': return False
    if smoothie4 != 'darkness': return False
    if animal2 != 'cat' or smoothie3 != 'desert': return False
    if child3 != 'billy': return False
    if name2 != 'carol' or animal3 != 'fish': return False
    if child2 != 'bella' or smoothie3 != 'butterscotch': return False
    if animal3 != 'dog': return False
    if child2 != 'bella' or child3 != 'timothy': return False
    
    return True

# Generate all permutations and check for valid arrangements
for names_perm in permutations(names):
    for children_perm in permutations(children):
        for animals_perm in permutations(animals):
            for smoothies_perm in permutations(smoothies):
                arrangement = {
                    'names': names_perm,
                    'children': children_perm,
                    'animals': animals_perm,
                    'smoothies': smoothies_perm
                }
                if is_valid(arrangement):
                    print(arrangement)
                    break