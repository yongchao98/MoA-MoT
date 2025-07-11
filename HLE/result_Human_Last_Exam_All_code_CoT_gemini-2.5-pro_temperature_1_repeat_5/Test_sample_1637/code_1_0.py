def calculate_hydrophobicity_score(molecule_name, fragments):
    """Calculates and prints a simple hydrophobicity score for a molecule."""
    
    # A simplified dictionary of hydrophobicity scores for molecular fragments.
    # Positive values are hydrophobic (poor water solubility).
    # Negative values are hydrophilic (good water solubility).
    score_map = {
        "Thioxanthenone Core": 5.0,
        "Chloro-hexyl Chain (-C6H12Cl)": 3.5,
        "Amide Group (-CONH-)": -1.0,
        "Short Ether Linker (-O-(CH2)2-O-(CH2)2-)": -1.5,
        "Acetate Linker (-O-CH2-)": 0.5,
        "PEG4 Chain (4x -CH2CH2O-)": -4.0, # Very hydrophilic
    }
    
    total_score = 0
    equation_parts = []
    
    print(f"--- Analysis for: {molecule_name} ---")
    for fragment in fragments:
        score = score_map.get(fragment, 0)
        total_score += score
        equation_parts.append(f"({fragment}: {score})")
        
    equation = " + ".join(equation_parts)
    print(f"Calculation: {equation}")
    print(f"Total Estimated Hydrophobicity Score: {total_score:.2f}\n")
    return total_score

def main():
    """Main function to run the analysis."""
    print("This script estimates the change in hydrophobicity for a chemical probe.\n")
    print("A lower score suggests better water solubility and less precipitation.\n")

    # Define fragments for the original probe
    original_probe_fragments = [
        "Thioxanthenone Core",
        "Acetate Linker (-O-CH2-)",
        "Amide Group (-CONH-)",
        "Short Ether Linker (-O-(CH2)2-O-(CH2)2-)",
        "Chloro-hexyl Chain (-C6H12Cl)"
    ]
    
    # Define fragments for the proposed new probe with a PEG chain
    # We replace the short ether and chloro-hexyl chain with a more soluble PEG chain.
    proposed_probe_fragments = [
        "Thioxanthenone Core",
        "Acetate Linker (-O-CH2-)",
        "Amide Group (-CONH-)",
        "PEG4 Chain (4x -CH2CH2O-)"
    ]

    original_score = calculate_hydrophobicity_score("Original Probe", original_probe_fragments)
    proposed_score = calculate_hydrophobicity_score("Proposed Probe with PEG", proposed_probe_fragments)
    
    if proposed_score < original_score:
        print("Conclusion: The proposed change significantly reduces the hydrophobicity score.")
        print("Therefore, replacing the hydrophobic alkyl chain with a hydrophilic PEG chain is very likely to solve the precipitation problem.")
    else:
        print("Conclusion: The proposed change does not appear to improve solubility based on this estimation.")

if __name__ == "__main__":
    main()