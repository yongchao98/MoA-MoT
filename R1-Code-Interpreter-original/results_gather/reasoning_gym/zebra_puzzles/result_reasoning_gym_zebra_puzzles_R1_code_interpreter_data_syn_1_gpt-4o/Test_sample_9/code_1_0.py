from itertools import permutations

names = ['bob', 'alice', 'arnold', 'carol']
colors = ['purple', 'white', 'brown', 'blue']
smoothies = ['dragonfruit', 'darkness', 'desert', 'butterscotch']
children = ['bella', 'billy', 'timothy', 'alice']

# We know some fixed positions:
# Bob loves purple
# Arnold is the mother of Billy and loves Dragonfruit
# The person who loves brown is in the second house

# Generate all permutations of the remaining attributes
for name_perm in permutations(names):
    for color_perm in permutations(colors):
        for smoothie_perm in permutations(smoothies):
            for child_perm in permutations(children):
                # Apply the known conditions
                if (name_perm[0] == 'bob' and color_perm[0] == 'purple' and
                    name_perm[2] == 'arnold' and child_perm[2] == 'billy' and smoothie_perm[2] == 'dragonfruit' and
                    color_perm[1] == 'brown' and
                    child_perm[0] == 'bella' and color_perm[0] == 'blue' and
                    child_perm[1] == 'alice' and name_perm[1] == 'alice' and
                    smoothie_perm[0] == 'butterscotch' and smoothie_perm[1] == 'desert' and
                    smoothie_perm[3] == 'darkness'):
                    # Check the left-right conditions
                    if (smoothie_perm[0] == 'butterscotch' and name_perm[1] == 'carol' and child_perm[3] == 'timothy' and
                        color_perm[0] == 'purple' and color_perm[1] == 'brown' and color_perm[2] == 'white' and color_perm[3] == 'blue'):
                        print(name_perm)
                        break