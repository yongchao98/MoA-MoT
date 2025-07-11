def find_product_iupac_name():
    """
    This script determines the IUPAC name for the product of the described reaction
    by applying chemical principles of Grignard reactions.
    """
    
    # 1. Define the parameters based on the reaction
    parent_molecule = "benzene"
    # The starting molecule has halogens at positions 1, 2, and 3.
    # The Grignard reagent provides a phenyl group as the substituent.
    substituent = "phenyl"
    # All halogens are replaced, so phenyl groups will be at these positions.
    positions = [1, 2, 3]

    # 2. Explain the chemical reasoning
    print("--- Chemical Reasoning ---")
    print("Reactant: 1,3-dibromo-2-iodobenzene")
    print("Reagent: Excess Phenyl Magnesium Bromide (a Grignard Reagent)")
    print("Reaction: The phenyl group from the Grignard reagent substitutes all halogen atoms on the benzene ring.")
    print("Product: A benzene molecule with phenyl groups at positions 1, 2, and 3.")
    print("-" * 26)

    # 3. Build the IUPAC name programmatically
    
    # Determine the correct numerical prefix (e.g., 'di', 'tri')
    count = len(positions)
    prefix_map = {2: "di", 3: "tri", 4: "tetra"}
    prefix = prefix_map.get(count, "")
    
    # Construct the position string, printing each number as required
    print("Constructing the name from the product structure:")
    position_string = ""
    for i, pos in enumerate(positions):
        print(f"Substituent position found at carbon: {pos}")
        position_string += str(pos)
        if i < len(positions) - 1:
            position_string += ","

    # Combine all parts into the final name
    final_iupac_name = f"{position_string}-{prefix}{substituent}{parent_molecule}"
    
    print("\n--- Final IUPAC Name ---")
    print(final_iupac_name)

# Run the function to get the answer
find_product_iupac_name()