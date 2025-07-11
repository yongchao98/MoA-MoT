def get_reaction_product():
    """
    Determines the IUPAC name of the product from the specified reaction.
    """
    
    # 1. Starting material and reaction analysis
    reaction = "Amide Enolate Ireland-Claisen Rearrangement"
    
    # 2. Deconstruction of the starting material for transformation mapping
    acyl_group = "propionyl"
    chiral_auxiliary = "N-((S)-1-phenylethyl)"
    allyl_group = "((S)-5-methylcyclopent-1-en-1-yl)methyl"

    # 3. Construction of the product's fragments based on the rearrangement mechanism
    parent_amide = "propanamide"
    
    # The spectator group remains unchanged
    product_N_substituent = chiral_auxiliary
    
    # The rearranged allyl group becomes a substituent on the C2 of the propanamide
    # Locants and stereochemistry are determined from reaction models
    rearranged_substituent_name = "((1S,4S)-4-methyl-2-methylenecyclopentan-1-yl)"
    
    # The stereochemistry of the new chiral center on the propanamide chain
    main_chain_stereochem = "(2R)"
    
    # 4. Assembling the final IUPAC name
    # Format: (Stereochem)-N-Substituent-C2_Substituent-Parent
    
    # We print the logic and numbers used to construct the name as requested.
    print("Logic for constructing the product's IUPAC name:")
    print(f"Reaction Type: {reaction}")
    print("The N-allyl group migrates to the C2 (alpha-carbon) of the propionyl group.")
    print("\n--- Name components ---")
    print(f"Parent Amide: {parent_amide}")
    print(f"Stereocenter at C2 of the main chain: {main_chain_stereochem}")
    print(f"Substituent on Nitrogen atom (unchanged): {chiral_auxiliary}")
    print(f"Substituent on Carbon 2 of the main chain: {rearranged_substituent_name}")
    print("Locants within the C2 substituent:")
    print("  - Point of attachment to main chain: Carbon 1")
    print("  - Exocyclic methylene (=CH2) group: Carbon 2")
    print("  - Methyl (CH3) group (stereocenter preserved): Carbon 4")
    
    # Combining the parts into the final name
    final_product_name = f"{main_chain_stereochem}-{product_N_substituent}-2-{rearranged_substituent_name}propanamide"

    print("\n--- Final Product IUPAC Name ---")
    print(final_product_name)


if __name__ == "__main__":
    get_reaction_product()
    # The final answer is wrapped in <<<>>> as requested.
    # Note: The stereochemistry of the new centers (2R, 1S) is the predicted major diastereomer
    # based on standard models for this highly stereoselective reaction.
    print("\n<<<(2R)-N-((S)-1-phenylethyl)-2-((1S,4S)-4-methyl-2-methylenecyclopentan-1-yl)propanamide>>>")