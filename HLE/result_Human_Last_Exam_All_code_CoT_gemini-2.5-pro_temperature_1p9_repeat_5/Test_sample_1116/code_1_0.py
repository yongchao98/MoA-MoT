import sys

# Disable writing __pycache__ files
sys.dont_write_bytecode = True

def analyze_bromination_product():
    """
    Analyzes a bromination reaction to identify the product based on H-NMR data.
    """

    # --- Step 1: Define Inputs ---
    starting_material_name = "2,8-bis(4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione"
    experimental_nmr_signals = 3

    # --- Step 2: Define Possible Products and their Properties ---
    # We analyze the three most likely species in the reaction mixture.
    possible_products = [
        {
            "name": "Starting Material",
            "description": "No reaction has occurred. Used as a baseline.",
            "symmetry": True,
            # Each of the 2 identical thiophenes has an H at positions 3 and 5.
            # Due to symmetry, there is one signal for the two H3 protons and one for the two H5 protons.
            "predicted_signals": 2,
            "reasoning": (
                "The starting material is symmetrical. It has two identical terminal thiophene rings. "
                "Each ring has two distinct aromatic protons (H3 and H5) which are coupled. "
                "Due to the molecule's symmetry, this results in a total of 2 signals in the aromatic region of the NMR spectrum."
            )
        },
        {
            "name": "Monobrominated Product",
            "full_name": "2-(4-(2-ethylhexyl)thiophen-2-yl)-8-(5-bromo-4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione",
            "description": "One of the two terminal thiophenes is brominated at its reactive C5 position.",
            "symmetry": False,
            # Unbrominated thiophene: H3 + H5 -> 2 signals
            # Brominated thiophene: H3' -> 1 signal
            "predicted_signals": 3,
            "reasoning": (
                "This molecule is asymmetrical. It has one unbrominated thiophene ring and one brominated ring.\n"
                "The unbrominated ring provides 2 signals (from its H3 and H5 protons).\n"
                "The brominated ring has lost its H5 proton and has only the H3' proton left, which gives 1 signal.\n"
                "The final equation for the total number of signals is: 2 + 1 = 3 signals."
            )
        },
        {
            "name": "Dibrominated Product",
            "description": "Both terminal thiophenes are brominated at their C5 positions.",
            "symmetry": True,
            # Both thiophenes are identical and have only an H3 proton left.
            # Due to symmetry, this results in only one signal for both H3 protons.
            "predicted_signals": 1,
            "reasoning": (
                "Symmetry is restored in this molecule. Both terminal thiophenes are now identically brominated at the C5 position. "
                "This leaves only one type of aromatic proton, H3. As both H3 protons are in identical chemical environments, "
                "they produce just 1 signal in the NMR spectrum."
            )
        }
    ]

    # --- Step 3: Analysis and Conclusion ---
    print("### Analysis of Bromination Reaction Product ###")
    print(f"\nStarting Material: {starting_material_name}")
    print(f"Experimental Observation: The new product shows {experimental_nmr_signals} peaks > 6.0 ppm in its ¹H-NMR spectrum.\n")
    print("--- Evaluating Possible Products ---")

    identified_product = None
    for product in possible_products:
        print(f"\nCandidate: {product['name']}")
        print(f"Predicted Aromatic Signals: {product['predicted_signals']}")
        print(f"Reasoning: {product['reasoning']}")
        
        if product['predicted_signals'] == experimental_nmr_signals:
            identified_product = product
            print(">>> MATCH: This prediction is consistent with the experimental data.")
        else:
            print(">>> MISMATCH: This prediction does not fit the experimental data.")

    print("\n--- Final Conclusion ---")
    if identified_product:
        print(f"The new spot isolated from the reaction is the '{identified_product['name']}'.")
        print("\nThis structure is the only one that explains the observation of exactly 3 peaks in the aromatic region of the ¹H-NMR spectrum.")
        print(f"\nFull Chemical Name: {identified_product.get('full_name', 'N/A')}")
    else:
        print("Could not definitively identify the product based on the provided NMR data.")

if __name__ == '__main__':
    analyze_bromination_product()