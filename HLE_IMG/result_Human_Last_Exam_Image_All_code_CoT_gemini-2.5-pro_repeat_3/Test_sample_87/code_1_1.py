def describe_products():
    """
    This function describes the chemical structures of products A, B, and C.
    """

    # Description of Product C
    print("Product C:")
    print("Molecular Formula: C11H16N2O3")
    print("Structure Description: Product C is the result of the N-acetylation of the secondary amine on the pyrrolidine ring of the starting material.")
    print("Systematic Name: N-acetyl-2-(1-pyrrolin-2-yl)pyrrolidine-2-carboxylic acid.")
    print("-" * 30)

    # Description of Product B
    print("Product B:")
    print("Molecular Formula: C12H14N2O3")
    print("Structure Description: Product B is formed in a two-step sequence. First, the starting material undergoes a Michael addition to methyl propiolate. The resulting intermediate then undergoes an intramolecular condensation between its carboxylic acid and methyl ester groups, eliminating methanol to form a bicyclic product containing a seven-membered cyclic anhydride ring.")
    print("-" * 30)

    # Description of Product A
    print("Product A:")
    print("Molecular Formula: C14H20N2O3")
    print("Structure Description: Product A is the major cycloadduct. First, the N-acetylated intermediate (Product C) undergoes decarboxylation (loss of CO2) to form a reactive azomethine ylide (a 1,3-dipole). This ylide is then trapped in a [3+2] dipolar cycloaddition reaction with a molecule of methyl propiolate, forming a tricyclic product.")
    print("-" * 30)

if __name__ == '__main__':
    describe_products()