def get_pummerer_product_name():
    """
    Determines the IUPAC name of the product from the reaction of methyl phenyl
    sulfoxide with triflic anhydride and trimethylsilyl cyanide.
    This script explains the reaction step-by-step.
    """

    # Step 1: Define reactants and reaction type
    reactant_sulfoxide = "Methyl phenyl sulfoxide"
    reactant_activator = "Triflic anhydride"
    reactant_nucleophile = "Trimethylsilyl cyanide"
    reaction_type = "Pummerer Reaction (cyanation variant)"

    print("--- Chemical Reaction Analysis ---")
    print(f"This is a {reaction_type}.")
    print(f"1. Substrate: {reactant_sulfoxide} (Ph-S(=O)-Me)")
    print(f"2. Activator: {reactant_activator} (Tf2O)")
    print(f"3. Nucleophile: {reactant_nucleophile} (TMSCN)")
    print("-" * 34 + "\n")

    # Step 2: Describe the reaction mechanism
    print("--- Reaction Mechanism ---")
    print("1. Activation: The sulfoxide oxygen attacks the highly electrophilic triflic anhydride to form a sulfoxonium salt intermediate, [Ph-S(OTf)-Me]+.")
    print("2. Thionium Ion Formation: A proton is lost from the acidic methyl group, and the triflate group (OTf-) leaves, forming the key Pummerer intermediate: a thionium ion, [Ph-S=CH2]+.")
    print("3. Nucleophilic Attack: Cyanide (from TMSCN) attacks the electrophilic carbon of the thionium ion. This forms the final carbon-cyanide bond.")
    print("-" * 26 + "\n")

    # Step 3: Identify the products and write the equation
    product_name = "2-(phenylthio)acetonitrile"
    product_structure = "C6H5-S-CH2-CN"
    byproduct1 = "Trimethylsilyl triflate (TMSOTf)"
    byproduct2 = "Triflic acid (TfOH)"

    print("--- Final Equation ---")
    # As requested, outputting each number (stoichiometric coefficient) in the final equation.
    stoichiometry = {
        reactant_sulfoxide: 1,
        reactant_activator: 1,
        reactant_nucleophile: 1,
        product_name: 1,
        byproduct1: 1,
        byproduct2: 1
    }

    equation_str = (
        f"{stoichiometry[reactant_sulfoxide]} {reactant_sulfoxide} + "
        f"{stoichiometry[reactant_activator]} {reactant_activator} + "
        f"{stoichiometry[reactant_nucleophile]} {reactant_nucleophile} -> \n"
        f"{stoichiometry[product_name]} {product_name} + "
        f"{stoichiometry[byproduct1]} {byproduct1} + "
        f"{stoichiometry[byproduct2]} {byproduct2}"
    )
    print(equation_str)
    print("-" * 22 + "\n")

    # Step 4: Final Answer
    print("The IUPAC name of the major organic product is:")
    print(product_name)

# Run the analysis
get_pummerer_product_name()
