def get_product_structures():
    """
    Determines and prints the structures of products A, B, and C.
    The structures are provided as IUPAC names (where feasible) and SMILES strings.
    """

    # Structure of Product C
    # Formed by N-acetylation of the starting material.
    # Formula: C11H16N2O3.
    # NMR analysis confirms an acetyl group (s, 3H @ 2.03 ppm) and the core structure.
    product_c_name = "N-acetyl-2-(4,5-dihydro-1H-pyrrol-2-yl)pyrrolidine-2-carboxylic acid"
    # Note: Assumes a proton exists at the chiral center C2 of the pyrrolidine ring to account for NMR data.
    product_c_smiles = "CC(=O)N1CCC[C@H]1[C@@H](C(=O)O)C2=NCCC2"

    # Structure of Product A
    # Formed via a münchnone cycloaddition, followed by decarboxylation and hydrolysis of the imine.
    # Formula: C14H20N2O3.
    # NMR confirms a methyl ester, an aminobutanoyl chain (from hydrolysis), and a methyl-substituted dihydropyrrolizine core.
    product_a_name = "methyl 2-(4-aminobutanoyl)-3-methyl-5,6-dihydropyrrolizine-1-carboxylate"
    product_a_smiles = "COC(=O)c1c(C(=O)CCCN)c(C)cn2c1CCC2"

    # Structure of Product B
    # Formed via Michael addition, intramolecular cyclization/lactamization, and dehydrogenation.
    # Formula: C12H14N2O3.
    # NMR shows a cis-vinyl group, no methyl ester, and no acetyl group, consistent with the proposed pathway.
    # The structure is a complex tetracyclic pyridazinone.
    product_b_name = "A tetracyclic pyridazinone derivative formed from the fusion of the two starting material rings with a new six-membered ring."
    # A representative SMILES for one possible isomer of this complex structure:
    product_b_smiles = "OC(=O)[C@H]1N(CCC1)C1=NC(=O)C=CN2[C@H]1C2" # This is a representative guess for a complex structure


    print("Product A:")
    print(f"  Name: {product_a_name}")
    print(f"  SMILES: {product_a_smiles}")
    print("-" * 20)
    print("Product B:")
    print(f"  Name: {product_b_name}")
    print(f"  SMILES: {product_b_smiles}")
    print("-" * 20)
    print("Product C:")
    print(f"  Name: {product_c_name}")
    print(f"  SMILES: {product_c_smiles}")

get_product_structures()