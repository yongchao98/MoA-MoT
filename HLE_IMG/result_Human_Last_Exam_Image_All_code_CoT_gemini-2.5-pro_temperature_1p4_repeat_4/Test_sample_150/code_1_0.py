def identify_reaction_product():
    """
    This function analyzes the provided multi-step chemical reaction and identifies
    the name of each intermediate and the final product.
    """
    
    steps = [
        {
            "step": 1,
            "reaction": "Benzene reacts with Propanoyl chloride in the presence of AlCl3 (Friedel-Crafts Acylation).",
            "product_name": "Intermediate-1",
            "product_structure": "Propiophenone"
        },
        {
            "step": 2,
            "reaction": "Propiophenone reacts with Br2/FeBr3 (Electrophilic Aromatic Bromination).",
            "product_name": "Intermediate-2",
            "product_structure": "3-Bromopropiophenone"
        },
        {
            "step": 3,
            "reaction": "3-Bromopropiophenone is reduced with H2/Pd (Catalytic Hydrogenation).",
            "product_name": "Intermediate-3",
            "product_structure": "3-Bromopropylbenzene"
        },
        {
            "step": 4,
            "reaction": "3-Bromopropylbenzene reacts with NBS and (PhCO2)2 (Radical Benzylic Bromination).",
            "product_name": "Final Product",
            "product_structure": "1-bromo-1-(3-bromophenyl)propane"
        }
    ]

    print("Analyzing the reaction sequence:")
    final_product = ""
    for s in steps:
        print(f"\nStep {s['step']}: {s['reaction']}")
        print(f"Resulting {s['product_name']}: {s['product_structure']}")
        if s['product_name'] == "Final Product":
            final_product = s['product_structure']
            
    print("\n-------------------------------------------")
    print(f"The systematic name of the final product is: {final_product}")
    print("-------------------------------------------")

identify_reaction_product()