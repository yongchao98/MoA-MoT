import textwrap

def solve_reaction_puzzle():
    """
    Identifies the product of a chemical reaction based on NMR data.
    """

    # 1. Define experimental observation
    experimental_nmr_peaks = 3

    # 2. Define potential compounds and their predicted NMR signals
    # based on molecular symmetry and number of aromatic protons.
    compounds = {
        "Starting Material": {
            "name": "2,8-bis(4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione",
            "reasoning": "Symmetrical, 2 equivalent outer thiophenes, each with 2 protons.",
            "predicted_nmr_peaks": 2
        },
        "Mono-brominated Product": {
            "name": "2-(5-bromo-4-(2-ethylhexyl)thiophen-2-yl)-8-(4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione",
            "reasoning": "Asymmetrical. One thiophene has 1 proton; the other has 2. Total of 3 non-equivalent protons.",
            "predicted_nmr_peaks": 3
        },
        "Di-brominated Product (Intended)": {
            "name": "2,8-bis(5-bromo-4-(2-ethylhexyl)thiophen-2-yl)-5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione",
            "reasoning": "Symmetrical. 2 equivalent outer thiophenes, each with 1 proton.",
            "predicted_nmr_peaks": 1
        }
    }

    # 3. Analyze and identify the correct product
    identified_product_name = "Unknown"
    
    print("--- Analysis of Potential Products vs. Experimental Data ---")
    print(f"Experimental ¹H-NMR data shows {experimental_nmr_peaks} peaks in the aromatic region (> 6.0 ppm).\n")

    for compound_type, properties in compounds.items():
        print(f"Candidate: {compound_type}")
        print(f"Predicted Aromatic Peaks: {properties['predicted_nmr_peaks']}")
        print(f"Reasoning: {properties['reasoning']}")
        if properties["predicted_nmr_peaks"] == experimental_nmr_peaks:
            identified_product_name = properties["name"]
            print("Status: MATCH FOUND\n")
        else:
            print("Status: No match\n")
            
    # 4. Print the final conclusion
    print("--- Conclusion ---")
    print("Based on the analysis, the new spot with 3 aromatic ¹H-NMR peaks corresponds to the mono-brominated product.")
    print("\nThe chemical name of the identified product is:")
    
    # Use textwrap to print the long name nicely
    wrapped_name = textwrap.fill(identified_product_name, width=70)
    print(wrapped_name)

# Run the analysis
solve_reaction_puzzle()
