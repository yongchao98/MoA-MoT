import re

def solve_reaction_sequence():
    """
    This script explains a 3-step chemical synthesis and identifies the final product.
    """

    # Define materials and products for clarity
    start_material = "(3S)-3-bromo-1-phenylbutane"
    product_A = "4-phenylbut-1-ene"
    product_B = "4-phenylbutan-1-ol"
    product_C = "1-bromo-4-phenylbutane"

    # Step 1: E2 Elimination
    print("--- Step 1: E2 Elimination ---")
    print(f"The starting material, {start_material}, is reacted with potassium tert-butoxide.")
    print("This strong, bulky base performs an E2 elimination, favoring the Hofmann product by abstracting a proton from the less hindered C4.")
    print(f"Product A is: {product_A}\n")

    # Step 2: Hydroboration-Oxidation
    print("--- Step 2: Hydroboration-Oxidation ---")
    print(f"Product A, {product_A}, is treated with borane followed by oxidative workup.")
    print("This is an anti-Markovnikov hydroboration-oxidation, adding an -OH group to the terminal carbon (C1) of the double bond.")
    print(f"Product B is: {product_B}\n")

    # Step 3: Bromination
    print("--- Step 3: Bromination ---")
    print(f"Product B, {product_B}, is treated with phosphorous tribromide (PBr3).")
    print("This reaction converts the primary alcohol into a primary alkyl bromide via an SN2 mechanism.")
    print(f"Product C is: {product_C}\n")

    # Final Product Analysis
    print("--- Final Product C Analysis ---")
    print(f"The final product is {product_C}.")
    print("IUPAC Name: 1-bromo-4-phenylbutane")
    print("Chirality: The starting material was chiral, but the stereocenter was removed in Step 1. No new stereocenters were formed. Thus, the final product is achiral.")

    # As requested, outputting the numbers from the final IUPAC name
    print("\nThe numbers in the final IUPAC name equation are:")
    numbers = re.findall(r'\d+', product_C)
    for num in numbers:
        print(num)

solve_reaction_sequence()
