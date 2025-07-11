import re

def get_wittig_product():
    """
    Determines the product of a Wittig reaction between pivalaldehyde and a specific phosphonium ylide.
    """
    # Step 1: Define the reactants by their functional components.
    # Pivalaldehyde: (CH3)3C-CHO
    aldehyde_alkyl_part = "(CH3)3C-CH" # The part that remains from the aldehyde
    aldehyde_carbonyl = "O"           # The oxygen that is replaced

    # Ylide: Ph3P=CH-CH2-(2-Cl-phenyl)
    ylide_alkylidene_part = "CH-CH2-(2-Cl-phenyl)" # The carbon part that replaces the oxygen
    ylide_phosphorus_part = "P(Ph)3"             # The phosphorus part that is eliminated

    # Step 2: Assemble the products based on the Wittig reaction mechanism.
    # The reaction forms an alkene and triphenylphosphine oxide.
    # Product: (CH3)3C-CH=CH-CH2-(2-Cl-phenyl)
    product_alkene_structure = f"{aldehyde_alkyl_part}={ylide_alkylidene_part}"

    # Byproduct: Ph3P=O
    byproduct_structure = "O=P(Ph)3"

    # Step 3: Determine the IUPAC name for the product.
    # Structure: (CH3)3C-CH=CH-CH2-(C6H4Cl)
    # The longest carbon chain containing the double bond is a 5-carbon chain (pentene).
    # C1 is part of the CH2 group attached to the ring. C2 and C3 form the double bond. C4 has two methyls.
    # This gives: 1-(2-chlorophenyl)-4,4-dimethylpent-2-ene.
    # Unstabilized ylides typically yield the (Z)-isomer.
    product_iupac_name = "1-(2-chlorophenyl)-4,4-dimethylpent-2-ene"
    byproduct_name = "triphenylphosphine oxide"

    # Step 4: Print the full reaction equation.
    # Use simpler, more common representations for the equation.
    pivalaldehyde_formula = "(CH3)3C-CHO"
    # The name "(2-(2-chlorophenyl)ethylidene)triphenyl-l5-phosphane" corresponds to the formula below
    ylide_formula = "Ph3P=CHCH2(2-Cl-C6H4)"
    product_formula = "(CH3)3C-CH=CH-CH2(2-Cl-C6H4)"
    byproduct_formula = "Ph3PO"

    print("The Wittig reaction involves the following transformation:\n")
    print(f"Aldehyde: {pivalaldehyde_formula}")
    print(f"Ylide: {ylide_formula}")
    print("-" * 40)
    print("The aldehyde's oxygen atom and the ylide's PPh3 group are eliminated and a C=C double bond is formed.\n")
    print("Products:")
    print(f"  Alkene: {product_formula}")
    print(f"  IUPAC Name: {product_iupac_name}")
    print(f"  Byproduct: {byproduct_formula} ({byproduct_name})\n")
    print("Final Equation:")
    # The problem asks to "output each number in the final equation!".
    # This likely means to show the full chemical formulas in the equation string.
    # The stoichiometric coefficients are all 1.
    final_equation = f"1 * [{pivalaldehyde_formula}] + 1 * [{ylide_formula}] -> 1 * [{product_formula}] + 1 * [{byproduct_formula}]"
    print(final_equation)


if __name__ == '__main__':
    get_wittig_product()