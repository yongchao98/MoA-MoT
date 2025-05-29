from itertools import permutations

# Define the possible attributes
names = ['arnold', 'bob', 'carol', 'alice']
nationalities = ['brit', 'chinese', 'dane', 'german']
colors = ['purple', 'white', 'blue', 'brown']
drinks = ['milk', 'water', 'boba tea', 'coffee']

# Generate all possible permutations of the attributes
name_perms = permutations(names)
nationality_perms = permutations(nationalities)
color_perms = permutations(colors)
drink_perms = permutations(drinks)

# Iterate over all possible combinations
for name_perm in name_perms:
    for nationality_perm in nationality_perms:
        for color_perm in color_perms:
            for drink_perm in drink_perms:
                # Apply the clues
                if (color_perm[nationality_perm.index('chinese')] == 'white' and
                    drink_perm[name_perm.index('arnold')] == 'coffee' and
                    nationality_perm.index('german') == nationality_perm.index('dane') - 1 and
                    color_perm[name_perm.index('carol')] == 'white' and
                    nationality_perm[nationality_perm.index('brit')] == color_perm.index('brown') and
                    color_perm[name_perm.index('arnold')] == 'purple' and
                    drink_perm[name_perm.index('bob')] == 'boba tea' and
                    name_perm.index('bob') == name_perm.index('alice') - 1 and
                    color_perm.index('purple') == drink_perm.index('milk') - 1):
                    
                    # If all conditions are met, print the arrangement
                    print("House 1:", name_perm[0], nationality_perm[0], color_perm[0], drink_perm[0])
                    print("House 2:", name_perm[1], nationality_perm[1], color_perm[1], drink_perm[1])
                    print("House 3:", name_perm[2], nationality_perm[2], color_perm[2], drink_perm[2])
                    print("House 4:", name_perm[3], nationality_perm[3], color_perm[3], drink_perm[3])