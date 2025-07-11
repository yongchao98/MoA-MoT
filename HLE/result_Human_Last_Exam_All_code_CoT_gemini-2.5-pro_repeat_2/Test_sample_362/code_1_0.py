def solve_wittig_reaction():
    """
    This function solves for the product of a Wittig reaction and prints the result.
    """
    # 1. Define the reactants
    aldehyde_name = "pivalaldehyde"
    aldehyde_carbonyl_part = "=O"
    aldehyde_alkyl_part = "(CH3)3C-CH"
    
    ylide_name = "(2-(2-chlorophenyl)ethylidene)triphenyl-l5-phosphane"
    ylide_phosphorus_part = "Ph3P="
    # The alkylidene part that replaces the carbonyl oxygen
    ylide_alkylidene_part = "CH-CH2-(C6H4-2-Cl)"

    # 2. Determine ylide type and resulting stereochemistry
    # The ylide has an alkyl group (not an electron-withdrawing group) attached to the P=C bond.
    # Therefore, it is a non-stabilized ylide.
    ylide_type = "non-stabilized"
    
    # Non-stabilized ylides predominantly form the Z-isomer.
    if ylide_type == "non-stabilized":
        stereochem_prefix = "(Z)-"
    else:
        # Stabilized ylides would give the E-isomer.
        stereochem_prefix = "(E)-"

    # 3. Construct the product and byproduct
    product_structure = f"{aldehyde_alkyl_part}={ylide_alkylidene_part}"
    product_iupac_name = "1-(2-chlorophenyl)-4,4-dimethylpent-2-ene"
    major_product_full_name = f"{stereochem_prefix}{product_iupac_name}"
    
    byproduct = "Ph3P=O (triphenylphosphine oxide)"

    # 4. Print the reaction equation and product information
    print("Wittig Reaction Analysis:")
    print("-" * 30)
    print(f"Aldehyde: {aldehyde_name} ({aldehyde_alkyl_part}{aldehyde_carbonyl_part})")
    print(f"Ylide: {ylide_name} ({ylide_phosphorus_part}{ylide_alkylidene_part})")
    print("-" * 30)
    
    print("\nReaction Equation:")
    print(f"{aldehyde_alkyl_part}{aldehyde_carbonyl_part}  +  {ylide_phosphorus_part}{ylide_alkylidene_part}  --->  {product_structure}  +  {byproduct}\n")
    
    print("Major Product Information:")
    print(f"Structure: {product_structure}")
    print(f"IUPAC Name: {major_product_full_name}")

    # As requested, output each number from the final product's IUPAC name
    product_numbers = {'phenyl_substituent_position': 1, 'chloro_position': 2, 'dimethyl_positions': [4, 4], 'alkene_position': 2}
    print("\nNumerical Locants in the Product Name:")
    print(f"Position of the (2-chlorophenyl) group: {product_numbers['phenyl_substituent_position']}")
    print(f"Position of the chloro group on the phenyl ring: {product_numbers['chloro_position']}")
    print(f"Position of the two methyl groups: {product_numbers['dimethyl_positions'][0]} and {product_numbers['dimethyl_positions'][1]}")
    print(f"Position of the double bond (pent-2-ene): {product_numbers['alkene_position']}")

# Execute the function
solve_wittig_reaction()
<<<(Z)-1-(2-chlorophenyl)-4,4-dimethylpent-2-ene>>>