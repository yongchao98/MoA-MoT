def identify_reaction_product():
    """
    This script analyzes a multi-step chemical reaction and identifies the final product.
    Each step is described, and the name of the compound at each stage is printed.
    """

    # --- Step 1 ---
    # Reaction: Friedel-Crafts Acylation
    # Benzene reacts with propanoyl chloride in the presence of AlCl3.
    # The acyl group (-COCH2CH3) attaches to the benzene ring.
    intermediate_1 = "1-phenylpropan-1-one"
    print("Step 1: The reaction between benzene and propanoyl chloride is a Friedel-Crafts Acylation.")
    print(f"          -> Intermediate-1 is: {intermediate_1}\n")

    # --- Step 2 ---
    # Reaction: Electrophilic Aromatic Bromination
    # 1-phenylpropan-1-one is reacted with Br2/FeBr3. The acyl group is a meta-director.
    # Bromine is added to the position meta to the acyl group.
    intermediate_2 = "1-(3-bromophenyl)propan-1-one"
    print("Step 2: The bromination of Intermediate-1 is an Electrophilic Aromatic Substitution.")
    print("          The acyl group directs the bromine to the meta position.")
    print(f"          -> Intermediate-2 is: {intermediate_2}\n")

    # --- Step 3 ---
    # Reaction: Catalytic Hydrogenation / Hydrogenolysis
    # Intermediate-2 is reduced with H2/Pd. This reduces the ketone to a methylene group (CH2)
    # and also removes the bromine atom (hydrogenolysis).
    intermediate_3 = "Propylbenzene"
    print("Step 3: The reaction with H2/Pd is a Catalytic Hydrogenation.")
    print("          Both the ketone and the C-Br bond are reduced.")
    print(f"          -> Intermediate-3 is: {intermediate_3}\n")

    # --- Step 4 ---
    # Reaction: Radical Benzylic Bromination
    # Propylbenzene is reacted with NBS and a radical initiator.
    # This specifically adds a bromine atom to the benzylic carbon (the carbon attached to the ring).
    final_product = "(1-bromopropyl)benzene"
    print("Step 4: The reaction with NBS is a Radical Benzylic Bromination.")
    print("          A bromine atom is added to the benzylic position.")
    print(f"          -> The Final Product is: {final_product}\n")


if __name__ == '__main__':
    identify_reaction_product()