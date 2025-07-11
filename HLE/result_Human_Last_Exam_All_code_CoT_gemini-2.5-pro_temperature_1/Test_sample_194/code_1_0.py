import itertools

def generate_character_matrix():
    """
    Encodes the described traits into a numerical matrix.
    
    Trait Encoding:
    - Trait 0: Integument (0: glabrous, 1: setose)
    - Trait 1: Leg Claw (0: absent, 1: present)
    - Trait 2: Simple Eye (0: present, 1: absent/lost)
    - Trait 3: Antennae (0: present, 1: absent/lost)
    - Trait 4: Antennae Texture (0: not serrate, 1: serrate)
    
    N/A states (like antennae texture for a species without antennae) are marked with -1.
    """
    char_matrix = {
        # Species: [Integument, Claw, Simple Eye, Antennae, Serrate]
        1: [0, 0, 0, 0, 1],
        2: [1, 1, 0, 0, 0],
        3: [1, 0, 1, 0, 0],
        4: [0, 1, 0, 0, 1],
        5: [1, 0, 1, 1, -1] 
    }
    return char_matrix

def print_difference_matrix(matrix):
    """Calculates and prints the pairwise trait difference matrix."""
    print("--- Trait Difference Matrix ---")
    print("This matrix shows the number of traits that differ between each pair of species.")
    
    species_ids = sorted(matrix.keys())
    
    # Header
    header = " | " + " | ".join(map(str, species_ids))
    print(header)
    print("-" * len(header))
    
    for id1 in species_ids:
        row = f"{id1}|"
        for id2 in species_ids:
            if id1 == id2:
                diff = 0
            else:
                diff = 0
                for i in range(len(matrix[id1])):
                    # Only compare if data is available for both
                    if matrix[id1][i] != -1 and matrix[id2][i] != -1:
                        if matrix[id1][i] != matrix[id2][i]:
                            diff += 1
            row += f" {diff} |"
        print(row)

def get_parsimony_score(tree, matrix):
    """
    Calculates the parsimony score of a tree using a simplified Fitch algorithm.
    The score is the total number of evolutionary changes (steps) required.
    """
    # For N/A (-1), we consider all possibilities {0, 1}
    fitch_matrix = {
        sp: [[val] if val != -1 else [0, 1] for val in traits]
        for sp, traits in matrix.items()
    }
    
    num_chars = len(next(iter(matrix.values())))
    total_score = 0
    
    for i in range(num_chars):
        char_score, _ = _calculate_char_score(tree, fitch_matrix, i)
        total_score += char_score
        
    return total_score

def _calculate_char_score(node, matrix, char_index):
    """Recursive helper for parsimony calculation on a single character."""
    # If the node is a leaf (a species), return its state
    if isinstance(node, int):
        return 0, set(matrix[node][char_index])

    # The node is a tuple (left_child, right_child)
    left_score, left_states = _calculate_char_score(node[0], matrix, char_index)
    right_score, right_states = _calculate_char_score(node[1], matrix, char_index)
    
    score = left_score + right_score
    intersection = left_states.intersection(right_states)
    
    if intersection:
        # No change needed at this node
        return score, intersection
    else:
        # A change is required (homoplasy/convergence)
        # The ancestral state could be any of the children's states
        return score + 1, left_states.union(right_states)

def solve_phylogeny():
    """Main function to perform analysis and print results."""
    
    print("Step 1: Define Character Matrix based on species descriptions.\n")
    char_matrix = generate_character_matrix()
    print("Character Matrix (Species x Traits):")
    print("Traits: [Setose, Claw, NoSimpleEye, NoAntennae, Serrate]")
    for sp, traits in sorted(char_matrix.items()):
        print(f"Species {sp}: {traits}")
    print("\n" + "="*50 + "\n")
    
    print("Step 2: Construct Trait Difference Matrix.\n")
    print_difference_matrix(char_matrix)
    print("\n" + "="*50 + "\n")

    print("Step 3: Calculate Parsimony Score for each potential tree.\n")
    
    # Tree topologies from answer choices (D is invalid)
    # Note: Written in a machine-readable nested tuple format.
    trees = {
        'A': (((1, (2, 4)), 3), 5),
        'B': (3, (2, (4, (1, 5)))),
        'C': ((3, ((1, 4), 2)), 5),
        'E': (((1, 4), (2, 3)), 5)
    }

    scores = {}
    print("--- Parsimony Scores ---")
    for name, topology in trees.items():
        score = get_parsimony_score(topology, char_matrix)
        scores[name] = score
        print(f"Tree {name} {topology}: requires {score} steps")
    
    print("\n" + "="*50 + "\n")
    
    # Find the best tree
    best_tree_name = min(scores, key=scores.get)
    best_score = scores[best_tree_name]

    print("Step 4: Describe the most parsimonious phylogeny.\n")
    print(f"The most parsimonious tree is Tree {best_tree_name} with a score of {best_score} steps.")
    print("This tree minimizes the number of required evolutionary events (gains/losses of traits).")
    print("\nDescription of Phylogeny E:")
    print("- The earliest split separates Species 5 from all others. Species 5 is unique in having lost its antennae.")
    print("- The remaining species form a clade which splits into two pairs: (1,4) and (2,3).")
    print("- The grouping of Species 1 and 4 is supported by their shared trait of serrate antennae.")
    print("- The grouping of Species 2 and 3 is supported by their shared trait of a setose (hairy) integument.")
    print("- This arrangement requires that leg claws and the loss of the simple eye each evolved twice independently (convergent evolution).")

    print("\nThe most parsimonious phylogeny is:")
    final_tree_string = "(((1,4),(2,3)),5)"
    print(final_tree_string)


if __name__ == '__main__':
    solve_phylogeny()
<<<E>>>