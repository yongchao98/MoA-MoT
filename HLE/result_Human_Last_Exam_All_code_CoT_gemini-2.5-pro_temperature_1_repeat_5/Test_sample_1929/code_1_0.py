def get_diels_alder_product():
    """
    This function describes the Diels-Alder reaction between butadiene and
    1,1-dichloro-2,2-difluoroethene and prints the reaction equation.
    """

    # --- Define Reactants ---
    reactant_1_name = "Butadiene"
    reactant_1_formula = "C4H6"

    reactant_2_name = "1,1-dichloro-2,2-difluoroethene"
    reactant_2_formula = "C2Cl2F2"

    # --- Determine Product ---
    # The reaction is a [4+2] cycloaddition which forms a cyclohexene ring.
    # Following IUPAC naming conventions for the resulting structure:
    product_name = "4,4-dichloro-5,5-difluorocyclohex-1-ene"
    product_formula = "C6H6Cl2F2"

    # --- Print the Reaction Equation ---
    print("The reaction between butadiene and 1,1-dichloro-2,2-difluoroethene is a Diels-Alder reaction.")
    print("The chemical equation for this reaction is:")
    
    # Print the equation, which includes all the required numbers in names and formulas.
    print(f"{reactant_1_name} ({reactant_1_formula}) + {reactant_2_name} ({reactant_2_formula}) -> {product_name} ({product_formula})")

# Execute the function to display the result
get_diels_alder_product()