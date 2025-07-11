def analyze_complex_stability():
    """
    Analyzes the stability of four iridium complexes to predict their relative lifetimes in LECs.
    """

    # Define the properties of each complex based on their structure
    complexes = {
        1: {"name": "Complex 1", "substituent": "2',4'-difluoro", "fluorine_count": 2},
        2: {"name": "Complex 2", "substituent": "4'-fluoro (para)", "fluorine_count": 1},
        3: {"name": "Complex 3", "substituent": "2'-fluoro (ortho)", "fluorine_count": 1},
        4: {"name": "Complex 4", "substituent": "none", "fluorine_count": 0}
    }

    print("Analysis of Complex Stability for LEC Lifetime Prediction:")
    print("-" * 60)
    print("Principle: The stability of these iridium complexes, and thus their device lifetime, is enhanced by fluorination of the phenylpyridine ligands. Fluorine strengthens the crucial Ir-C bond.")
    print("\nRanking the complexes by expected stability (from most to least stable):")

    print("\n1. Complex 1 (Most Stable):")
    print("   - Has two fluorine atoms per ligand. The strong, combined electron-withdrawing effect provides maximum stabilization to the Ir-C bond, leading to the longest expected lifetime.")

    print("\n2. Complex 2 (Stable):")
    print("   - Has one fluorine at the para-position. This provides significant electronic stabilization without introducing steric strain, leading to a long lifetime.")

    print("\n3. Complex 3 (Less Stable):")
    print("   - Has one fluorine at the ortho-position. While it provides electronic stabilization, its proximity to the metal center can introduce steric strain/repulsion. This potential instability makes it less stable than Complex 2.")

    print("\n4. Complex 4 (Least Stable):")
    print("   - Has no fluorine atoms. It lacks the electronic stabilization of the other complexes, making its Ir-C bond the weakest. It is expected to have the shortest lifetime.")

    print("-" * 60)
    print("\nConclusion:")
    print("The question asks for the complexes expected to show SHORTER lifetimes.")
    print("Based on the stability ranking (1 > 2 > 3 > 4), the two least stable complexes are 3 and 4.")

    shorter_lifetime_complexes = [3, 4]
    print(f"\nFinal Answer: The complexes with shorter lifetimes are [{shorter_lifetime_complexes[0]}, {shorter_lifetime_complexes[1]}]")


analyze_complex_stability()