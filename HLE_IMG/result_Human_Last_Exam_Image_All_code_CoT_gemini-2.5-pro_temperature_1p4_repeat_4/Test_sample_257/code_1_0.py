def solve_nmr_problem():
    """
    This script determines the splitting pattern and integration of the most deshielded proton
    in Compound 1 by analyzing the chemical reaction and molecular structure.
    """

    # --- Step 1: Analyze the structure and identify the proton of interest ---
    # The reaction is a sulfonation, which adds -SO3H groups to the outer rings
    # of the Pr-DAOTA cation. This doesn't change the protons on the central ring.
    # The most deshielded proton in this aromatic system is the one on the central
    # ring, situated between the two outer rings (in the "bay region").
    
    # --- Step 2: Determine the integration of the signal for this proton ---
    # The molecule retains its C2 symmetry after the reaction.
    # There is only one proton in this specific chemical environment.
    integration = 1

    # --- Step 3: Determine the splitting pattern for this proton's signal ---
    # The splitting pattern follows the n+1 rule, where 'n' is the number of
    # equivalent neighboring protons.
    # This proton is on the central ring and has two neighbors on the same ring.
    # Due to the molecule's C2 symmetry, these two neighbors are chemically equivalent.
    num_equivalent_neighbors = 2
    
    # The number of peaks in the multiplet is n + 1.
    num_peaks = num_equivalent_neighbors + 1
    
    if num_peaks == 1:
        splitting_pattern = "singlet"
    elif num_peaks == 2:
        splitting_pattern = "doublet"
    elif num_peaks == 3:
        splitting_pattern = "triplet"
    elif num_peaks == 4:
        splitting_pattern = "quartet"
    else:
        splitting_pattern = f"multiplet with {num_peaks} peaks"

    # --- Step 4: Print the final conclusion ---
    print("Conclusion for the highest deshielded proton peak in Compound 1:")
    print(f"Integration: The peak corresponds to {integration} proton.")
    print(f"Splitting Pattern: It is a {splitting_pattern}, as it is coupled to {num_equivalent_neighbors} equivalent neighboring protons, resulting in {num_peaks} peaks.")

solve_nmr_problem()