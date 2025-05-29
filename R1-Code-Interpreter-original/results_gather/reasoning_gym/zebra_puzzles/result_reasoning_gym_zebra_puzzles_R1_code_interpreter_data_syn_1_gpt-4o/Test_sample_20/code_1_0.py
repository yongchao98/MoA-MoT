from itertools import permutations

# Possible names, nationalities, lunches, and smoothies
names = ['alice', 'bob', 'carol', 'arnold']
nationalities = ['chinese', 'german', 'dane', 'brit']
lunches = ['pizza', 'grilled cheese', 'soup', 'stir fry']
smoothies = ['desert', 'butterscotch', 'darkness', 'dragonfruit']

# Fixed positions based on clues
fixed_positions = {
    'bob': 0,  # Bob is in the first house
    'carol': 3,  # Carol is in the fourth house
    'dane': 0,  # Dane is in the first house
    'chinese': names.index('alice'),  # Alice is Chinese
}

# Check all permutations
for name_perm in permutations(names):
    if name_perm[0] != 'bob' or name_perm[3] != 'carol':
        continue
    for nat_perm in permutations(nationalities):
        if nat_perm[0] != 'dane' or nat_perm[names.index('alice')] != 'chinese':
            continue
        for lunch_perm in permutations(lunches):
            if lunch_perm[nat_perm.index('german')] != 'pizza':
                continue
            for smoothie_perm in permutations(smoothies):
                if smoothie_perm[lunch_perm.index('pizza')] != 'desert':
                    continue
                if smoothie_perm[lunch_perm.index('stir fry')] != 'darkness':
                    continue
                if smoothie_perm[nat_perm.index('brit')] != 'dragonfruit':
                    continue
                if nat_perm.index('brit') + 1 != lunch_perm.index('grilled cheese'):
                    continue
                # If all conditions are satisfied, print the result
                print(name_perm, nat_perm, lunch_perm, smoothie_perm)