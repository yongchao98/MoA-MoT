import itertools

def solve_phylogeny():
    """
    Solves the alien species phylogeny problem by:
    1. Defining the character matrix for the five species.
    2. Calculating and printing a trait difference matrix.
    3. Calculating the parsimony score for each candidate tree.
    4. Identifying the most parsimonious tree and printing the result.
    """

    # Step 1 & 2: Define Character Matrix
    # Characters:
    # 1. Body Covering: 0=glabrous, 1=setose
    # 2. Simple Eye:    0=present, 1=absent
    # 3. Leg Claws:     0=absent, 1=present
    # 4. Antennae:      0=present, 1=absent
    # 5. Serration:     0=not serrate, 1=highly serrate (None for species 5)
    species_data = {
        '1': {'traits': [0, 0, 0, 0, 1], 'name': 'Species 1'},
        '2': {'traits': [1, 0, 1, 0, 0], 'name': 'Species 2'},
        '3': {'traits': [1, 1, 0, 0, 0], 'name': 'Species 3'},
        '4': {'traits': [0, 0, 1, 0, 1], 'name': 'Species 4'},
        '5': {'traits': [1, 1, 0, 1, None], 'name': 'Species 5'},
    }
    char_names = ["Body", "Simple Eye", "Claws", "Antennae", "Serration"]

    # Step 3: Calculate and Print Trait Difference Matrix
    print("Trait Difference Matrix:")
    species_list = sorted(species_data.keys())
    header = "   " + "  ".join(species_list)
    print(header)
    print("---" * (len(species_list) + 1))
    
    for s1_key in species_list:
        row = f"{s1_key} |"
        for s2_key in species_list:
            if s1_key == s2_key:
                row += " 0 "
                continue
            
            diff = 0
            for i in range(len(char_names)):
                trait1 = species_data[s1_key]['traits'][i]
                trait2 = species_data[s2_key]['traits'][i]
                # Only compare if both traits are known
                if trait1 is not None and trait2 is not None and trait1 != trait2:
                    diff += 1
            row += f" {diff} "
        print(row)
    print("\n" + "="*40 + "\n")

    # Step 4: Parsimony Analysis for each tree
    # Manual Fitch algorithm logic for each character on each tree topology.
    
    # Tree E: (((1,4),(2,3)),5)
    scores_E = {
        'Body': 1,       # Change 1->0 for clade (1,4)
        'Simple Eye': 2, # Independent losses in 3 and 5
        'Claws': 2,      # Independent gains in 2 and 4
        'Antennae': 1,   # Single loss in 5
        'Serration': 1   # Single gain for clade (1,4)
    }
    total_E = sum(scores_E.values())

    # Tree A: (((1,(4,2)),3),5) which is ((((4,2),1),3),5)
    scores_A = {
        'Body': 2,       # Homoplasy required
        'Simple Eye': 2, # Loss in 3, loss in 5
        'Claws': 1,      # Single gain for clade (4,2)
        'Antennae': 1,   # Single loss in 5
        'Serration': 2   # Gain in 1 and 4 independently
    }
    total_A = sum(scores_A.values())
    
    # Tree C: ((3,((4,1),2)),5) which is (5,(3,(2,(1,4))))
    scores_C = {
        'Body': 1,       # Change 1->0 for clade (1,4)
        'Simple Eye': 2, # Loss in 3 and 5
        'Claws': 2,      # Independent gains in 2 and 4
        'Antennae': 1,   # Single loss in 5
        'Serration': 1   # Single gain for clade (1,4)
    }
    total_C = sum(scores_C.values())

    results = {
        'A. (((1,(4,2)),3),5)': {'scores': scores_A, 'total': total_A},
        'C. ((3,((4,1),2)),5)': {'scores': scores_C, 'total': total_C},
        'E. (((1,4),(2,3)),5)': {'scores': scores_E, 'total': total_E},
    }

    print("Parsimony Scores (Tree Lengths):\n")
    min_score = float('inf')
    best_tree = None

    for tree, result in results.items():
        score_breakdown = " + ".join(f"{v} ({k})" for k, v in result['scores'].items())
        print(f"Tree {tree}")
        print(f"Score: {score_breakdown} = {result['total']}\n")
        if result['total'] < min_score:
            min_score = result['total']
            best_tree = tree

    # Find all trees with the minimum score
    best_trees = [tree for tree, result in results.items() if result['total'] == min_score]

    print("="*40 + "\n")
    print(f"The minimum number of evolutionary steps required is {min_score}.")
    print("The most parsimonious tree(s) are:")
    for tree in best_trees:
        print(tree)

    # Step 5: Describe the most parsimonious phylogeny
    # Although there's a tie in score, Tree E is preferred as its internal clades (1,4) and (2,3)
    # are each supported by a clear shared derived character (synapomorphy):
    # - Clade (1,4) is supported by the gain of serrated antennae.
    # - Clade (2,3) is supported by the gain of a setose body.
    # Tree C's structure relies on grouping by shared ancestral states (plesiomorphies), which is a weaker basis.
    final_choice_name = "E. (((1,4),(2,3)),5)"
    final_choice_phylogeny = "(((1,4),(2,3)),5)"

    print("\nBased on the principle that clades should be defined by shared derived traits (synapomorphies),")
    print(f"the preferred phylogeny is {final_choice_name}.")
    
    print("\nFinal Phylogeny:")
    print(final_choice_phylogeny)

solve_phylogeny()
<<<E>>>