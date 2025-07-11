def solve_chemical_reaction():
    """
    This script tracks the transformations in a multi-step organic synthesis
    to identify the name of the final product.
    """
    # Step 1: Friedel-Crafts Acylation
    reactant_1 = "Benzene"
    reagent_1 = "Propanoyl chloride / AlCl3"
    intermediate_1 = "Propiophenone (1-phenylpropan-1-one)"
    print("Step 1: Friedel-Crafts Acylation")
    print(f"   Reactant: {reactant_1}")
    print(f"   Reagent: {reagent_1}")
    print(f"   Product (Intermediate-1): {intermediate_1}")
    print("-" * 50)

    # Step 2: Electrophilic Aromatic Bromination
    reactant_2 = intermediate_1
    reagent_2 = "Br2 / FeBr3"
    intermediate_2 = "1-(3-bromophenyl)propan-1-one"
    print("Step 2: Electrophilic Bromination")
    print(f"   Reactant: {reactant_2}")
    print(f"   Reagent: {reagent_2}")
    print(f"   Description: The acyl group is a meta-director.")
    print(f"   Product (Intermediate-2): {intermediate_2}")
    print("-" * 50)

    # Step 3: Catalytic Hydrogenation
    reactant_3 = intermediate_2
    reagent_3 = "H2 / Pd"
    intermediate_3 = "Propylbenzene"
    print("Step 3: Catalytic Reduction")
    print(f"   Reactant: {reactant_3}")
    print(f"   Reagent: {reagent_3}")
    print(f"   Description: The ketone is reduced to an alkane and the aryl bromide is removed (hydrogenolysis).")
    print(f"   Product (Intermediate-3): {intermediate_3}")
    print("-" * 50)

    # Step 4: Free-Radical Bromination
    reactant_4 = intermediate_3
    reagent_4 = "NBS / (PhCO2)2, CCl4"
    final_product = "(1-Bromopropyl)benzene"
    print("Step 4: Free-Radical Benzylic Bromination")
    print(f"   Reactant: {reactant_4}")
    print(f"   Reagent: {reagent_4}")
    print(f"   Description: NBS selectively brominates the most stable radical position, which is the benzylic position.")
    print(f"   Product (Final Product): {final_product}")
    print("-" * 50)

if __name__ == "__main__":
    solve_chemical_reaction()