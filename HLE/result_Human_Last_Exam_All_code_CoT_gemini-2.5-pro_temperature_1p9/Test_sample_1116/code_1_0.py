def analyze_bromination_reaction():
    """
    Analyzes the reaction based on NMR data to propose a product structure.
    """

    # --- Define Molecule Properties ---
    # The starting material is symmetrical and has three distinct aromatic proton environments.
    sm_name = "2,8-bis(4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione"
    sm_signals = {
        "core_protons_H1_H7": 2,  # 2H, one signal due to symmetry
        "thiophene_beta_protons_H3": 2,  # 2H, one signal due to symmetry
        "thiophene_alpha_protons_H5": 2   # 2H, one signal due to symmetry
    }
    num_sm_signals = len(sm_signals)

    # --- The Puzzle ---
    # The key observation is that the product has the same number of aromatic signals as the starting material, but is a new compound.
    num_product_signals_observed = 3

    # --- Proposing the Product Structure ---
    # Standard alpha-bromination of the thiophenes would result in a product with 2 signals, which contradicts the observation.
    # Hypothesis: Radical bromination occurred on the 'benzylic' position of the 2-ethylhexyl side chains.
    # This reaction preserves the molecule's symmetry and all its aromatic protons.
    product_name = "2,8-bis(4-(1-bromo-2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione"
    product_signals = {
        "core_protons_H1_H7": 2,        # Unchanged, 2H, one signal
        "thiophene_beta_protons_H3": 2, # Unchanged, 2H, one signal
        "thiophene_alpha_protons_H5": 2  # Unchanged, 2H, one signal
    }
    num_proposed_product_signals = len(product_signals)


    # --- Output the Reasoning and Final Answer ---
    print("--- Analysis of the Chemical Puzzle ---")
    print(f"Starting Material: {sm_name}")
    print(f"Number of aromatic signals observed in starting material: {num_sm_signals}")
    print("This corresponds to 3 pairs of equivalent protons.\n")

    print(f"Observation: A new compound was isolated with {num_product_signals_observed} aromatic NMR signals.")
    print("This rules out the expected product from thiophene alpha-bromination, which would only have 2 signals.\n")

    print("--- Proposed Identity of the New Spot ---")
    print("The product is likely formed by radical bromination on the side chains, not the aromatic rings.")
    print(f"Proposed Product: {product_name}")
    print("This reaction creates a new molecule while preserving the symmetry and the three distinct sets of aromatic protons.\n")

    print("--- Final Equation of Proton Signals ---")
    print("This section shows the number of protons corresponding to each signal set in the final proposed product.")
    
    # "output each number in the final equation"
    # The final 'equation' is the description of the product's NMR signals.
    signal_count = 1
    for signal_name, proton_count in product_signals.items():
        print(f"Signal {signal_count}: {signal_name}, integrates to {proton_count} protons.")
        signal_count += 1
    
    total_protons = sum(product_signals.values())
    print(f"\nTotal Aromatic Protons = {total_protons}")
    print(f"Total Number of Aromatic Signals = {num_proposed_product_signals}")


if __name__ == '__main__':
    analyze_bromination_reaction()
<<<2,8-bis(4-(1-bromo-2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione>>>