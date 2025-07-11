import textwrap

class Molecule:
    """A class to represent a molecule and its predicted NMR properties."""
    def __init__(self, name, description, aromatic_proton_signals):
        self.name = name
        self.description = description
        self.aromatic_proton_signals = aromatic_proton_signals

    def print_details(self):
        """Prints the details of the molecule."""
        print(f"Name: {self.name}")
        print(f"Description: {self.description}")
        print(f"Predicted Aromatic 1H-NMR Signals (> 6.0 ppm): {self.aromatic_proton_signals}\n")

def solve_bromination_puzzle():
    """
    Solves the chemical puzzle by analyzing the starting material, reaction,
    and potential products to identify the compound matching the NMR data.
    """
    print("--- Step 1: Analysis of the Starting Material ---")
    # The starting material (SM) is symmetric. The two `(4-(2-ethylhexyl)thiophen-2-yl)` groups are identical.
    # Each outer thiophene has two protons in the aromatic region:
    # 1. H-5 (alpha to sulfur): Very reactive to electrophilic substitution.
    # 2. H-3 (beta to sulfur).
    # Due to the molecule's symmetry, there are only two distinct aromatic proton environments.
    sm = Molecule(
        name="2,8-bis(4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione",
        description="Symmetric molecule. The two outer thiophene rings are equivalent.",
        aromatic_proton_signals=2
    )
    sm.print_details()

    print("--- Step 2: Analysis of Potential Products ---")
    print("The reaction is an electrophilic bromination with NBS, which targets the most electron-rich positions.")
    print("The most reactive sites are the C5 positions (alpha-protons) of the outer thiophene rings.\n")

    # Potential Product 1: Monobrominated product
    # Brominating one of the two C5 positions makes the molecule asymmetric.
    # The three remaining aromatic protons are now all in unique chemical environments:
    # 1. H-5 on the unreacted thiophene.
    # 2. H-3 on the unreacted thiophene.
    # 3. H-3 on the brominated thiophene.
    # This results in 3 expected NMR signals.
    product_mono_bromo = Molecule(
        name="2-(5-bromo-4-(2-ethylhexyl)thiophen-2-yl)-8-(4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione",
        description="Asymmetric product from bromination at one C5 position.",
        aromatic_proton_signals=3
    )
    print("Candidate 1: Monobrominated Product")
    product_mono_bromo.print_details()

    # Potential Product 2: Dibrominated product
    # Brominating both C5 positions results in a symmetric molecule again.
    # Only the two H-3 protons remain, and they are chemically equivalent.
    # This results in 1 expected NMR signal.
    product_di_bromo = Molecule(
        name="2,8-bis(5-bromo-4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione",
        description="Symmetric product from bromination at both C5 positions. This is the expected final product.",
        aromatic_proton_signals=1
    )
    print("Candidate 2: Dibrominated Product")
    product_di_bromo.print_details()

    print("--- Step 3: Conclusion based on Experimental Data ---")
    experimental_peaks = 3
    print(f"The experimental observation is that the isolated new product has {experimental_peaks} peaks > 6.0 ppm in its H-NMR spectrum.")
    print("The sluggish reaction (initially no product, then formation after 1 hr) suggests the reaction may be incomplete.\n")

    # Compare candidates with experimental data
    candidates = [product_mono_bromo, product_di_bromo]
    final_product = None
    for candidate in candidates:
        if candidate.aromatic_proton_signals == experimental_peaks:
            final_product = candidate
            break

    if final_product:
        print("--- Final Answer ---")
        print("The product that matches the experimental data is the monobrominated compound, as it is the only candidate expected to have 3 aromatic H-NMR signals.")
        
        print("\nFinal Chemical Equation:")
        print("\nReactant:")
        print(textwrap.fill(sm.name, 80))
        print("\n  + NBS (in excess, 1 hr) --->")
        print("\nProduct:")
        print(textwrap.fill(final_product.name, 80))
        
        # This is the final answer for the <<<>>> block
        global final_answer_name
        final_answer_name = final_product.name
    else:
        print("Error: Could not identify the product.")
        final_answer_name = "Unknown"

if __name__ == '__main__':
    final_answer_name = ""
    solve_bromination_puzzle()
    # The final answer is printed here in the required format, though not part of the script's visible output
    # print(f"\n<<<{final_answer_name}>>>")
