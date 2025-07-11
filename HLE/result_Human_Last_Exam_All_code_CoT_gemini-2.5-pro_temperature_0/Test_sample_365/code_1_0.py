def get_product_formula():
    """
    Calculates and prints the molecular formula of the product.

    The reaction is an anionic oxy-Cope rearrangement, which is an intramolecular
    rearrangement. Therefore, the product has the same molecular formula as the
    starting material.

    The starting material is:
    (1S,2R,4S)-2-((S)-4-((tert-butyldimethylsilyl)oxy)cyclopent-1-en-1-yl)-7,7-dimethoxybicyclo[2.2.1]hept-5-en-2-ol

    Calculation of the molecular formula:
    - Core (7,7-dimethoxybicyclo[2.2.1]hept-5-en-2-ol): C9H14O3
    - Substituent ((S)-4-((tert-butyldimethylsilyl)oxy)cyclopent-1-en-1-yl): C11H21OSi
    - Combining them removes one H from each, resulting in:
      C = 9 + 11 = 20
      H = 13 + 21 = 34
      O = 3 + 1 = 4
      Si = 1
    """
    
    # The "equation" is: Starting Material -> Product. The formula represents the product.
    product_formula = {
        "Carbon (C)": 20,
        "Hydrogen (H)": 34,
        "Oxygen (O)": 4,
        "Silicon (Si)": 1
    }

    print("The reaction is an anionic oxy-Cope rearrangement.")
    print("The product is a bicyclic ketone with the same molecular formula as the starting material.")
    print("\nThe molecular formula of the product is C20H34O4Si.")
    print("\nThe elemental composition (the numbers in the final formula) is:")
    for element, count in product_formula.items():
        print(f"- {element}: {count}")

get_product_formula()