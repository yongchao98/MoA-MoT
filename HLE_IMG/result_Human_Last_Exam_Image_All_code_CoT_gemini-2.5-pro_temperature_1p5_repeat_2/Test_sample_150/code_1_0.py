def identify_product():
    """
    This function tracks the transformations in a multi-step organic synthesis
    to identify the name of the final product.
    """

    print("--- Reaction Sequence Analysis ---")

    # Step 1: Friedel-Crafts Acylation
    step1_desc = "Step 1 involves the Friedel-Crafts acylation of benzene with propanoyl chloride and AlCl3."
    intermediate1_name = "Propiophenone (1-phenylpropan-1-one)"
    print(step1_desc)
    print("Intermediate-1 is: {}\n".format(intermediate1_name))

    # Step 2: Electrophilic Aromatic Bromination
    step2_desc = "Step 2 is the bromination of propiophenone using Br2/FeBr3. The acyl group is a meta-director."
    intermediate2_name = "1-(3-bromophenyl)propan-1-one"
    print(step2_desc)
    print("Intermediate-2 is: {}\n".format(intermediate2_name))

    # Step 3: Catalytic Hydrogenation
    step3_desc = "Step 3 is the catalytic hydrogenation with H2/Pd. This reduces the ketone to an alkane and also removes the bromine from the aromatic ring (hydrogenolysis)."
    intermediate3_name = "Propylbenzene"
    print(step3_desc)
    print("Intermediate-3 is: {}\n".format(intermediate3_name))

    # Step 4: Radical Benzylic Bromination
    step4_desc = "Step 4 is the benzylic bromination of propylbenzene using NBS and a radical initiator. Bromine adds to the carbon attached to the benzene ring."
    final_product_name = "1-bromo-1-phenylpropane"
    print(step4_desc)
    print("The Final Product is: {}".format(final_product_name))

identify_product()