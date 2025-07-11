def solve_coordination_chemistry():
    """
    This script determines the atoms coordinated to a metal center based on
    the reactants in a coordination chemistry problem.
    """

    # Step 1: Define the properties of the ligand based on its name.
    # Ligand: 1,3-Di[3-(2’-pyridyl)pyrazol-1-ylmethyl]-2,4,6-triethylbenzene
    # 'Di' means two arms. Each arm has a pyridyl (1 N donor) and a pyrazolyl (1 N donor).
    ligand_donor_atoms = ['N', 'N', 'N', 'N']
    
    # Step 2: Define the properties of the metal salt.
    # Metal Salt: ZnBr2
    # Metal center is Zn(II). Anions are two Br- ions.
    # Zn(II) is a d10 metal and commonly forms 6-coordinate complexes.
    metal_preferred_coordination_number = 6
    anionic_ligands = ['Br', 'Br']

    # Step 3: Predict the final coordination sphere.
    # The strong, chelating ligand will coordinate first.
    coordinated_atoms = list(ligand_donor_atoms)
    
    # Calculate the number of remaining sites to fill to reach the preferred coordination number.
    sites_to_fill = metal_preferred_coordination_number - len(coordinated_atoms)
    
    # Fill the remaining sites with the next best available ligands, which are the bromide anions.
    # The solvent (methanol) is a weak ligand and is unlikely to coordinate.
    if sites_to_fill > 0:
        num_anions_to_add = min(sites_to_fill, len(anionic_ligands))
        for i in range(num_anions_to_add):
            coordinated_atoms.append(anionic_ligands[i])

    # Sort the list for consistent ordering (e.g., alphabetically).
    coordinated_atoms.sort()

    print("Based on the chemical principles:")
    print(f"1. The ligand is tetradentate, providing {len(ligand_donor_atoms)} Nitrogen donors.")
    print(f"2. The metal Zn(II) prefers a coordination number of {metal_preferred_coordination_number}.")
    print(f"3. The remaining {sites_to_fill} sites are filled by the {len(anionic_ligands)} Bromide anions.")
    print("\nTherefore, the atoms coordinated to the Zn center are:")
    
    # To match the answer format (B, B, N, N, N, N), we can print in that order.
    final_output_list = []
    for atom_type in ['Br', 'N']:
        count = coordinated_atoms.count(atom_type)
        for _ in range(count):
            final_output_list.append(atom_type)
            
    print(', '.join(final_output_list))

solve_coordination_chemistry()