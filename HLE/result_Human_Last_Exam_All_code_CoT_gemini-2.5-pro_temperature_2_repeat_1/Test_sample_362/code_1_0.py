def solve_wittig_reaction():
    """
    Determines and explains the product of a specific Wittig reaction.
    """
    # 1. Define Reactants
    aldehyde = {
        "name": "pivalaldehyde",
        "structure": "(CH3)3C-CHO"
    }
    ylide = {
        "name": "(2-(2-chlorophenyl)ethylidene)triphenyl-l5-phosphane",
        "structure": "Ph3P=CH-CH2-(o-Cl-C6H4)"
    }

    # 2. Define Products
    # The reaction replaces the =O of the aldehyde with the =C< part of the ylide.
    # Structure: (CH3)3C-CH=CH-CH2-(o-Cl-C6H4)
    # The bulky tert-butyl and substituted benzyl groups favor the E (trans) isomer.
    product_name = "(E)-1-(2-chlorophenyl)-4,4-dimethylpent-2-ene"
    side_product_name = "triphenylphosphine oxide"

    # 3. Construct and Print the Full Reaction Equation
    # The final equation includes names with locant numbers.
    equation = (
        f"{aldehyde['name']} + {ylide['name']} -->\n"
        f"{product_name} + {side_product_name}"
    )

    print("--- Wittig Reaction Analysis ---")
    print("\nThe full reaction equation is:\n")
    print(equation)

    # 4. Explain the naming of the main product
    print("\n--- Product Name Breakdown ---")
    print(f"The principal organic product is: {product_name}\n")
    print("The numbers in the name have the following meaning:")
    print(" - 'pent-2-ene': Indicates a 5-carbon chain with a double bond between carbon number 2 and carbon number 3.")
    print(" - '1-(2-chlorophenyl)': Specifies that a '(2-chlorophenyl)' group is attached to carbon number 1 of the chain.")
    print(" - '4,4-dimethyl': Indicates that two methyl groups are attached to carbon number 4.")
    print(" - '(E)-': Specifies the stereochemistry of the double bond is 'trans', which is the major product due to steric hindrance.")

# Execute the function to print the solution
solve_wittig_reaction()