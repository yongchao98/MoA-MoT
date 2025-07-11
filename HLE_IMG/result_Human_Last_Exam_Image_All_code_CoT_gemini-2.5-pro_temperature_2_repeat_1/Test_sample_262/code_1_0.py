def solve_lec_lifetime():
    """
    Identifies which Iridium complexes are expected to have shorter lifetimes in LECs.
    """
    # Step 1: Define the complexes based on the substitution pattern of their
    # cyclometalating (C^N) ligands. We focus on the positions of Fluorine atoms
    # on the phenyl ring, relative to the Ir-C bond (position 1').
    complexes = {
        1: {'name': 'Complex 1', 'f_positions': [2, 4], 'ligand_name': '2-(2,4-difluorophenyl)pyridine'},
        2: {'name': 'Complex 2', 'f_positions': [], 'ligand_name': '2-phenylpyridine'},
        3: {'name': 'Complex 3', 'f_positions': [2], 'ligand_name': '2-(2-fluorophenyl)pyridine'},
        4: {'name': 'Complex 4', 'f_positions': [4], 'ligand_name': '2-(4-fluorophenyl)pyridine'}
    }

    # Step 2: State the guiding chemical principle for shorter lifetimes in this context.
    print("--- Analysis of Emitter Lifetime in Light-Emitting Electrochemical Cells (LECs) ---")
    print("Principle: Iridium(III) complexes with fluorine atoms at the 2'-position (ortho to the Ir-C bond)")
    print("are known to be less stable under the reductive conditions present in an operating LEC.")
    print("This leads to C-F bond cleavage, a major degradation pathway that results in shorter device lifetimes.\n")

    # Step 3: Iterate through the complexes and identify those with the destabilizing feature.
    shorter_lifetime_complexes = []
    print("Screening complexes for ortho-fluorine (2'-F) substituents:")
    for complex_id, properties in complexes.items():
        # Check if 2 is in the list of fluorine positions
        has_ortho_fluorine = 2 in properties['f_positions']
        
        if has_ortho_fluorine:
            status = "HAS ortho-F, expecting SHORTER lifetime."
            shorter_lifetime_complexes.append(complex_id)
        else:
            status = "Does NOT have ortho-F, expecting longer lifetime."
            
        print(f" - Complex {complex_id} (Ligand: {properties['ligand_name']}): {status}")

    # Step 4: Present the conclusion.
    shorter_lifetime_complexes.sort()
    
    print("\n--- Conclusion ---")
    print("The complexes expected to show shorter lifetimes are those identified above.")
    
    # Per instruction, output the numbers in the final result.
    equation_str = " + ".join(map(str, shorter_lifetime_complexes))
    print(f"The equation for the numbers of the complexes with shorter lifetimes is: {equation_str}")
    
    final_set = sorted(shorter_lifetime_complexes)
    print(f"Final list of complexes with shorter lifetimes: {final_set}")


solve_lec_lifetime()