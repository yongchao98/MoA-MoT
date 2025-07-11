#
# Chemical Reasoning to Determine the Coordination Sphere of a Zinc Complex
#

def solve_coordination_chemistry():
    """
    This function deduces the atoms coordinated to a metal center based on
    the reactants and principles of coordination chemistry.
    """

    # 1. Define the reactants
    # The ligand is 1,3-Di[3-(2’-pyridyl)pyrazol-1-ylmethyl]-2,4,6-triethylbenzene.
    # Each arm has one pyridyl nitrogen and one pyrazolyl nitrogen that can coordinate.
    # Since there are two arms ('Di-'), the ligand is tetradentate.
    ligand_donor_atoms = ['N', 'N', 'N', 'N']

    # The metal salt is Zinc Bromide (ZnBr2).
    metal_center = 'Zn'
    salt_anions = ['Br', 'Br']

    # The solvent is methanol, which has a potential oxygen donor atom.
    solvent_donor_atom = 'O'

    # 2. Determine the coordination
    # The tetradentate ligand will bind to the zinc center.
    coordinated_atoms = list(ligand_donor_atoms)
    print(f"The tetradentate ligand provides {len(coordinated_atoms)} nitrogen atoms: {', '.join(coordinated_atoms)}")

    # Zinc(II) commonly has a coordination number of 4, 5, or 6.
    # With a tetradentate ligand, it will likely form a 5- or 6-coordinate complex.
    # The remaining coordination sites will be filled by the best available ligands.
    # We compare the anionic bromide (Br-) with the neutral solvent methanol (O donor).
    # Anionic ligands like bromide are generally preferred over neutral solvent molecules.
    coordinated_atoms.extend(salt_anions)
    print(f"The two bromide anions from ZnBr2 also coordinate.")

    # 3. Final Result
    # The final coordination sphere is formed by the 4 N atoms and the 2 Br atoms.
    # This results in a stable 6-coordinate complex.
    coordinated_atoms.sort()
    
    print("\nFinal coordinated atoms around the Zn center are:")
    # Using a loop to print each atom as requested in the prompt
    for atom in coordinated_atoms:
        print(atom)
        
    # The final answer corresponds to the choice 'Br, Br, N, N, N, N'.

solve_coordination_chemistry()