def solve_wittig_reaction():
    """
    Determines the product of a Wittig reaction between pivalaldehyde and
    (2-(2-chlorophenyl)ethylidene)triphenyl-l5-phosphane, and prints the detailed derivation.
    """
    # 1. Define Reactants and their fragments
    aldehyde_name = "pivalaldehyde"
    aldehyde_structure_fragment = "(CH3)3C-CH="
    
    ylide_name = "(2-(2-chlorophenyl)ethylidene)triphenyl-l5-phosphane"
    ylide_structure_fragment = "=CH-CH2-(2-chlorophenyl)"

    print("--- The Wittig Reaction ---")
    print(f"This reaction couples an aldehyde with a phosphorus ylide to form an alkene.")
    print("\nStep 1: Identify the reacting fragments.")
    print(f"Aldehyde: {aldehyde_name}")
    print(f"This provides the fragment: {aldehyde_structure_fragment}")
    print(f"\nYlide: {ylide_name}")
    print(f"This provides the fragment: {ylide_structure_fragment}")

    print("\nStep 2: Combine fragments to form the product backbone.")
    print("The C=O and C=P bonds are replaced by a new C=C double bond.")
    product_backbone_structure = "(CH3)3C-CH=CH-CH2-(2-chlorophenyl)"
    print(f"Combined structure: {product_backbone_structure}")
    
    print("\nStep 3: Determine the systematic IUPAC name for the product.")
    parent_chain_name = "pent-2-ene"
    locant_double_bond = 2
    substituent1 = "1-(2-chlorophenyl)"
    locant_phenyl_group = 1
    locant_chloro_group = 2
    substituent2 = "4,4-dimethyl"
    locant_methyl_groups = "4, 4"
    
    print(f"The longest carbon chain containing the double bond is a pentene: 5 carbons.")
    print(f"Numbering starts from the carbon attached to the phenyl ring to give the double bond the lowest number ({locant_double_bond}).")
    print(f"Base Name: {parent_chain_name}")
    print(f"Substituents:")
    print(f"- A (2-chlorophenyl) group at position {locant_phenyl_group}")
    print(f"- Two methyl groups at position(s) {locant_methyl_groups}")
    
    iupac_base_name = "1-(2-chlorophenyl)-4,4-dimethylpent-2-ene"
    print(f"\nPutting it together alphabetically gives the base name: {iupac_base_name}")

    print("\nStep 4: Determine the stereochemistry of the double bond.")
    print("The ylide is non-stabilized because the carbon alpha to the phosphorus is not adjacent to an electron-withdrawing group like a carbonyl.")
    print("Wittig reactions with non-stabilized ylides predominantly yield the Z-alkene (cis-isomer).")
    stereo_prefix = "(Z)"

    final_product_name = f"{stereo_prefix}-{iupac_base_name}"
    
    # Z-isomer SMILES for 1-(2-chlorophenyl)-4,4-dimethylpent-2-ene
    # The directional bonds / and \ indicate the relative stereochemistry.
    # CC(C)(C)/C=C\C... means the tert-butyl and the next carbon are on the same side ("cis").
    final_product_smiles = "CC(C)(C)/C=C\\Cc1c(Cl)cccc1"
    
    print("\n--- Final Product ---")
    print(f"Name: {final_product_name}")
    print(f"SMILES Representation: {final_product_smiles}")
    
    print("\n--- Numbers in Final IUPAC Name ---")
    print(f"Position of the chlorophenyl group: {locant_phenyl_group}")
    print(f"Position of the chloro on the phenyl ring: {locant_chloro_group}")
    print(f"Position of the double bond: {locant_double_bond}")
    print(f"Positions of the methyl groups: {locant_methyl_groups}")

if __name__ == '__main__':
    solve_wittig_reaction()