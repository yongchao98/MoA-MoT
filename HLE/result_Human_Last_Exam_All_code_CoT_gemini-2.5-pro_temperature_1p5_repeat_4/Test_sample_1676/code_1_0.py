def identify_compound_3():
    """
    This script analyzes a three-step chemical synthesis to identify the final product, Compound 3.
    It prints the step-by-step transformation, including the chemical formulas and atom counts for the final reaction.
    """
    
    print("--- Analysis of the Reaction Sequence ---")

    # Step 1: Epoxidation of Terpinolene
    print("\nStep 1: Terpinolene (C10H16) reacts with m-CPBA.")
    print("Result: Selective epoxidation of the endocyclic double bond.")
    print("Compound 1 is Terpinolene 1,2-oxide (C10H16O).")

    # Step 2: Conversion to Episulfide
    print("\nStep 2: Compound 1 (epoxide) reacts with N,N-dimethyl thioformamide and acid.")
    print("Result: The epoxide is converted to an episulfide (thiirane).")
    print("Compound 2 is Terpinolene 1,2-episulfide (C10H16S).")

    # Step 3: Reduction of Episulfide
    print("\nStep 3: Compound 2 (episulfide) is reduced with LiAlH4.")
    print("Result: Reductive opening of the episulfide ring to form a thiol.")
    
    # Final reaction equation and atom counts, as requested.
    reactant_formula = "C10H16S"
    product_formula = "C10H18S"
    final_product_name = "1-methyl-4-(propan-2-ylidene)cyclohexane-1-thiol"
    
    print(f"\nThe final chemical equation is: {reactant_formula} + 2[H] -> {product_formula}")
    
    print(f"\nAtom count for the reactant ({reactant_formula}):")
    print("Carbon (C): 10")
    print("Hydrogen (H): 16")
    print("Sulfur (S): 1")
    
    print(f"\nAtom count for the product ({product_formula}):")
    print("Carbon (C): 10")
    print("Hydrogen (H): 18")
    print("Sulfur (S): 1")
    
    print("\n--- Final Result ---")
    print(f"The chemical structure of Compound 3 is: {final_product_name}")

if __name__ == '__main__':
    identify_compound_3()