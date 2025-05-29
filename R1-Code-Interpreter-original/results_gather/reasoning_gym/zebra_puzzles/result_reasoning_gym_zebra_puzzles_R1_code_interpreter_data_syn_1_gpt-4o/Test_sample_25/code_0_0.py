from itertools import permutations

# Define the possible options
names = ['alice', 'arnold', 'carol', 'bob']
phones = ['google pixel 6', 'oneplus 9', 'huawei p50', 'samsung galaxy s21']
children = ['alice', 'billy', 'timothy', 'bella']
flowers = ['lilies', 'daffodils', 'iris', 'carnations']

# Iterate over all permutations of the arrangements
for name_perm in permutations(names):
    for phone_perm in permutations(phones):
        for child_perm in permutations(children):
            for flower_perm in permutations(flowers):
                # Apply the clues
                if (name_perm[0] == 'alice' and phone_perm[0] == 'google pixel 6' and
                    child_perm[3] == 'timothy' and phone_perm[3] == 'samsung galaxy s21' and
                    child_perm[1] == 'billy' and
                    flower_perm[1] == 'lilies' and phone_perm[1] == 'oneplus 9' and
                    flower_perm[2] == 'daffodils' and
                    name_perm[1] == 'arnold' and name_perm[2] == 'carol' and
                    flower_perm[0] == 'iris' and
                    child_perm[0] == 'alice' and phone_perm[2] == 'huawei p50'):
                    # If all conditions are met, print the name of the person in house 1
                    print(name_perm[0])