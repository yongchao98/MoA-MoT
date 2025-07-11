import re

def identify_elimination_product():
    """
    Analyzes the reaction of (1S,2R)-1-bromo-2-methylcyclohexane with potassium tert-butoxide
    and identifies the major product.
    """

    reactant_name = "(1S,2R)-1-bromo-2-methylcyclohexane"
    reagent = "potassium tert-butoxide (a strong, bulky base)"

    print("--- Reaction Analysis ---")
    print(f"The reaction involves treating {reactant_name} with {reagent}.")
    print("This is a classic E2 (bimolecular elimination) reaction.")
    print("\n")

    print("--- Stereochemistry of E2 on Cyclohexane ---")
    print("A critical requirement for the E2 mechanism on a cyclohexane ring is that the leaving group (Br) and a beta-proton (H on an adjacent carbon) must be in a trans-diaxial orientation (180 degrees apart).")
    print("\n")

    print("--- Substrate Conformation ---")
    print(f"In {reactant_name}, the '(1S,2R)' configuration means the bromo and methyl groups are cis to each other.")
    print("Of the two possible chair conformations, only the one where the Bromine atom is in an axial position can undergo the E2 reaction.")
    print("In this reactive conformer, the bulkier methyl group is in the more stable equatorial position.")
    print("\n")

    print("--- Regioselectivity: Zaitsev vs. Hofmann ---")
    print("From the reactive conformer (axial Br), elimination can occur in two ways:")
    print("1. Removing a proton from C2 (the carbon with the methyl group) yields 1-methylcyclohexene (Zaitsev product).")
    print("2. Removing a proton from C6 yields 3-methylcyclohexene (Hofmann product).")
    print("\n")
    print("Potassium tert-butoxide is a sterically bulky base. It has difficulty accessing the hindered proton on C2. It will preferentially attack the more accessible proton on C6.")
    print("Therefore, the reaction favors the Hofmann rule, forming the less substituted alkene as the major product.")
    print("\n")

    print("--- Final Product ---")
    product_name = "3-methylcyclohexene"
    print(f"The major product of the reaction is: {product_name}")
    print("\n")
    
    print("--- Numbers from the Chemical Equation ---")
    print(f"Equation: {reactant_name} --> {product_name}")
    
    reactant_numbers = re.findall(r'\d+', reactant_name)
    product_numbers = re.findall(r'\d+', product_name)

    print("Reactant numbers:", " ".join(reactant_numbers))
    print("Product numbers:", " ".join(product_numbers))


identify_elimination_product()
<<<3-methylcyclohexene>>>