import textwrap

def identify_bromination_product():
    """
    This script analyzes a chemical reaction to identify an unknown product based on NMR data.
    """
    # Define the molecular components for clarity
    core = "5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione"
    substituent = "4-(2-ethylhexyl)thiophen-2-yl"
    brominated_substituent = "5-bromo-4-(2-ethylhexyl)thiophen-2-yl"

    # --- Step 1: Analyze possible products and their expected NMR signals ---
    print("Step 1: Analyzing ¹H-NMR Signals for Potential Products")
    print("The key to solving this is to count the number of unique aromatic protons for each possible product.")
    print("-" * 70)

    # Starting Material
    print("A) Starting Material:")
    print(f"   Name: 2,8-bis({substituent})-{core}")
    print("   Symmetry: Symmetric")
    print("   Aromatic Protons: Two equivalent protons at C3 positions and two equivalent protons at C5 positions.")
    print("   Expected number of aromatic peaks > 6.0 ppm: 2")
    print("")

    # Di-bromo Product (fully reacted)
    print("B) Di-bromo Product (Symmetric):")
    print(f"   Name: 2,8-bis({brominated_substituent})-{core}")
    print("   Symmetry: Symmetric")
    print("   Aromatic Protons: Both C5 protons are replaced by Br. Only the two equivalent C3 protons remain.")
    print("   Expected number of aromatic peaks > 6.0 ppm: 1")
    print("")

    # Mono-bromo Product (intermediate)
    print("C) Mono-bromo Product (Asymmetric):")
    product_name = f"2-({substituent})-8-({brominated_substituent})-{core}"
    print(f"   Name: {textwrap.fill(product_name, width=67, initial_indent=' ' * 3, subsequent_indent=' ' * 9)}")
    print("   Symmetry: Asymmetric")
    print("   Aromatic Protons: One thiophene has H at C3 and C5. The other has only H at C3.")
    print("   Expected number of aromatic peaks > 6.0 ppm: 3")
    print("-" * 70)

    # --- Step 2: Compare with observation and conclude ---
    print("Step 2: Conclusion")
    print("The experimental data reports that the new isolated spot has THREE peaks larger than 6.0 ppm.")
    print("This observation perfectly matches the predicted spectrum for the asymmetric MONO-BROMO PRODUCT (C).")
    print("\nThe product is the result of adding one bromine atom to one of the outer thiophene rings.")
    print("-" * 70)

    # --- Step 3: Display the final deduced reaction equation ---
    print("Step 3: The Deduced Reaction Equation")
    print("\nReactant (Expected 2 aromatic peaks):")
    print(f"  2,8-bis({substituent})-{core}")
    print("\n   + ~1 eq. NBS reacts\n   --------------------->")
    print("\nProduct (Observed 3 aromatic peaks):")
    print(f"  {textwrap.fill(product_name, width=67, initial_indent='', subsequent_indent='  ')}")

identify_bromination_product()

final_answer = f"2-({ '4-(2-ethylhexyl)thiophen-2-yl' })-8-({ '5-bromo-4-(2-ethylhexyl)thiophen-2-yl' })-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione"
print(f"\n\n<<<What is this new spot?>>>\n<<<{final_answer}>>>")