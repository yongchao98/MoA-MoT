def solve_chemistry_problem():
    """
    This function analyzes the stability of four Iridium complexes
    to determine which ones are expected to have shorter lifetimes in LECs.
    """
    # Define the ligands for each complex. We use a simplified notation.
    # 'ppy' = 2-phenylpyridine (unsubstituted, less stable)
    # 'fppy' = fluorinated 2-phenylpyridine (more stable)
    complexes_data = {
        1: {'name': 'Complex 1', 'ligands': ['fppy', 'fppy']},
        2: {'name': 'Complex 2', 'ligands': ['fppy', 'ppy']},
        3: {'name': 'Complex 3', 'ligands': ['fppy', 'fppy']},
        4: {'name': 'Complex 4', 'ligands': ['fppy', 'ppy']}
    }

    # The presence of an unsubstituted 'ppy' ligand is the key factor for shorter lifetime.
    unstable_ligand_type = 'ppy'
    shorter_lifetime_complexes = []

    print("Analyzing complexes for the presence of the less stable ligand (unsubstituted 2-phenylpyridine)...")

    for complex_id, data in complexes_data.items():
        if unstable_ligand_type in data['ligands']:
            shorter_lifetime_complexes.append(complex_id)
            print(f"- {data['name']} (ID: {complex_id}) contains the unstable ligand and is expected to have a shorter lifetime.")
        else:
            print(f"- {data['name']} (ID: {complex_id}) does not contain the unstable ligand and is expected to have a longer lifetime.")

    shorter_lifetime_complexes.sort()
    
    print("\n--- Conclusion ---")
    print("The degradation of Iridium(III) emitters is often initiated at the site of the least stable ligand.")
    print("Fluorination of the C^N ligand strengthens the Ir-C bond, increasing stability and device lifetime.")
    print("Complexes with an unsubstituted 2-phenylpyridine (ppy) ligand are less stable.")
    print("\nThe complexes expected to show shorter lifetimes are:")
    # Print the final result clearly as requested
    print(f"[{shorter_lifetime_complexes[0]}, {shorter_lifetime_complexes[1]}]")

solve_chemistry_problem()