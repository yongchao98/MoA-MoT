import pandas as pd

def identify_starting_material():
    """
    Identifies a starting material by logically filtering candidates
    based on NMR spectroscopic data.
    """

    # Step 1: Define possible starting materials and their expected NMR signatures
    # after a two-step disubstitution.
    candidates = {
        '2,4-dichloropyrimidine': {'protons': 'H-5, H-6', 'coupling': 'coupled (doublets)'},
        '2,5-dichloropyrimidine': {'protons': 'H-4, H-6', 'coupling': 'coupled (doublets)'},
        '4,5-dichloropyrimidine': {'protons': 'H-2, H-6', 'coupling': 'uncoupled (singlets)'},
        '4,6-dichloropyrimidine': {'protons': 'H-2, H-5', 'coupling': 'uncoupled (singlets)'},
    }

    print("Step 1: Analyzing 1H NMR data to determine the substitution pattern.")
    # The 1H NMR spectrum shows two 1H singlets at 8.24 ppm and 8.11 ppm.
    # This indicates the product has two uncoupled protons on its core aromatic ring.
    h1_observation = "uncoupled (singlets)"
    print(f"Observation: Product's 1H NMR shows two aromatic protons with a coupling pattern of '{h1_observation}'.\n")
    
    # Filter candidates based on the 1H NMR coupling pattern.
    possible_candidates = {name: data for name, data in candidates.items() if data['coupling'] == h1_observation}
    
    print("Filtering candidates based on 1H NMR data...")
    print("The following candidates match the uncoupled proton pattern:")
    for name, data in possible_candidates.items():
        print(f"- {name} (Expected protons: {data['protons']})")
    print("-" * 20)

    # Step 2: Analyze 13C NMR data to further differentiate the remaining candidates.
    print("\nStep 2: Analyzing 13C NMR data to refine the structure.")
    c13_shifts_core = [156.89, 154.96, 152.80, 102.23] # Excluding substituent peaks
    print(f"Core 13C NMR signals (in ppm): {c13_shifts_core}")

    # Theoretical analysis of carbon environments for the remaining candidates.
    # For a 4,6-disubstituted pyrimidine, we expect C-H signals for C-2 and C-5. C-5 is
    # electronically shielded, appearing around 100-110 ppm.
    # For a 4,5-disubstituted pyrimidine, we expect C-H signals for C-2 and C-6. Both
    # are relatively deshielded and a highly shielded signal like 102.23 ppm is not expected.
    c13_analysis = {
        '4,6-dichloropyrimidine': "Product has C-H carbons at positions 2 and 5. C-5 is characteristically shielded (~100-110 ppm).",
        '4,5-dichloropyrimidine': "Product has C-H carbons at positions 2 and 6. Both are relatively deshielded (>130 ppm)."
    }
    
    print("\nEvaluating remaining candidates against 13C NMR data:")
    final_candidate = None
    for name in possible_candidates:
        print(f"Candidate: {name}")
        print(f"Analysis: {c13_analysis[name]}")
        # Check if the candidate's expected C13 signature matches the observation.
        # The key observation is the presence of a highly shielded carbon at 102.23 ppm.
        if "shielded (~100-110 ppm)" in c13_analysis[name] and 100 < 102.23 < 110:
            print("Result: Consistent with observed signal at 102.23 ppm.")
            final_candidate = name
        else:
            print("Result: Inconsistent with observed 13C NMR data.")
        print("-" * 20)

    # Final conclusion
    print("\nConclusion:")
    if final_candidate:
        print(f"The evidence from both 1H and 13C NMR strongly supports that Compound A is {final_candidate}.")
    else:
        print("Could not determine the starting material from the provided data.")

# Run the analysis
if __name__ == "__main__":
    identify_starting_material()
