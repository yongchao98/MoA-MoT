def print_structure_descriptions():
    """
    Prints the descriptions and SMILES representations of the three chemical structures A, B, and C.
    """
    
    # Description for Product A (C14H20N2O3)
    desc_A = (
        "Product A is the [3+2] cycloadduct formed between methyl propiolate and an N-acetylated azomethine ylide.\n"
        "The reaction proceeds in steps:\n"
        "1. First, an N-acetylated azomethine ylide is generated via thermal decarboxylation of the N-acetylated starting material (Product C).\n"
        "2. This ylide then undergoes a [3+2] cycloaddition with methyl propiolate to form Product A.\n"
        "The resulting structure is a tricyclic pyrrolizidine core, fused to the N-acetylated imidazoline ring, with a methyl ester substituent on the newly formed ring.\n"
        "SMILES representation for visualization: CC(=O)N1CCN=C(C1)C23CCCN3C=C(C2)C(=O)OC"
    )

    # Description for Product B (C12H14N2O3)
    desc_B = (
        "Product B is a complex tetracyclic dihydropyridone derivative formed via a multi-step cascade reaction.\n"
        "The reaction proceeds in steps:\n"
        "1. A non-acetylated azomethine ylide is generated from the starting material via decarboxylation.\n"
        "2. This ylide undergoes a [3+2] cycloaddition with methyl propiolate.\n"
        "3. The resulting adduct undergoes an intramolecular Michael addition.\n"
        "4. The tetracyclic intermediate is then oxidized (loses 4 hydrogens).\n"
        "5. Finally, the amidine functional group is hydrated to an amide, yielding Product B.\n"
        "SMILES representation for visualization: COC(=O)C1=CN2C(=O)[C@H]3N(CCC3)[C@H]2C=C1"

    )

    # Description for Product C (C11H16N2O3)
    desc_C = (
        "Product C is the N-acetylated derivative of the starting material.\n"
        "The reaction involves the acetylation of the secondary amine within the 4,5-dihydro-1H-imidazole ring by acetic anhydride.\n"
        "The pyrrolidine and carboxylic acid portions of the molecule remain unchanged.\n"
        "SMILES representation for visualization: CC(=O)N1CCN=C(C1)C2(C(=O)O)CCCN2"
    )

    print("--- Structure of Product A (C14H20N2O3) ---")
    print(desc_A)
    print("\n" + "="*60 + "\n")
    print("--- Structure of Product B (C12H14N2O3) ---")
    print(desc_B)
    print("\n" + "="*60 + "\n")
    print("--- Structure of Product C (C11H16N2O3) ---")
    print(desc_C)

print_structure_descriptions()