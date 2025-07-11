import sys

def generate_liquid_crystal_candidates():
    """
    Generates candidate liquid crystal molecules based on the structure:
    Alkyl-Chain -- Benzene-Ring -- Cyano-Group
    This script varies the alkyl chain length 'n' to explore potential
    candidates for room temperature liquid crystals, as per the design strategy.
    """
    # The polar head group is cyano (-CN)
    polar_group_smiles = "(C#N)"
    
    # The core is a 1,4-substituted benzene ring
    # The SMILES c1cc(X)ccc1Y represents para-substitution
    # We will build the full SMILES string inside the loop
    
    print("Generating candidate molecules for the series: 4-alkyl-1-cyanobenzene\n")
    print("-" * 50)
    print(f"{'Chain Length (n)':<20}{'Molecule Name':<25}{'SMILES String'}")
    print("-" * 50)

    # The tuning strategy suggests varying the alkyl chain length (n).
    # We will explore a range of n from 3 to 8, which includes the
    # recommended starting point of n=5.
    min_n = 3
    max_n = 8

    for n in range(min_n, max_n + 1):
        # Generate the alkyl chain part of the name
        # n=3: Propyl, n=4: Butyl, n=5: Pentyl, etc.
        prefixes = {3: 'Propyl', 4: 'Butyl', 5: 'Pentyl', 6: 'Hexyl', 7: 'Heptyl', 8: 'Octyl'}
        name_alkyl_part = prefixes.get(n, f'C{n}-alkyl')
        molecule_name = f"4-{name_alkyl_part}benzonitrile"

        # Generate the alkyl chain SMILES string (e.g., n=5 is CCCCC)
        alkyl_chain_smiles = "C" * n

        # Combine the parts to form the full SMILES string for the molecule
        # R-c1ccc(X)cc1 --> [alkyl]c1ccc([CN])cc1
        smiles_string = f"{alkyl_chain_smiles}c1ccc{polar_group_smiles}cc1"

        # Print the results in a formatted table
        print(f"{n:<20}{molecule_name:<25}{smiles_string}")

if __name__ == '__main__':
    generate_liquid_crystal_candidates()
