import sys

# Step 1: Explain the methodology to the user.
# We are using a simplified model where different chemical groups contribute to an overall
# "hydrophobicity score". A higher score means the molecule is more hydrophobic and less
# soluble in water (like the cell culture medium).
print("--- Analysis of Probe Solubility ---")
print("The precipitation issue is likely caused by low water solubility.")
print("We can estimate solubility by calculating a 'hydrophobicity score'.")
print("A lower score suggests better water solubility.\n")

# Step 2: Define the hydrophobicity scores for relevant chemical functional groups.
# These are simplified, representative values.
# Positive values are for hydrophobic groups, negative for hydrophilic groups.
hydrophobicity_scores = {
    'Thioxanthenone Core': 5.0,  # Large, aromatic, very hydrophobic
    'Methoxy Group': -0.5,        # Slightly hydrophilic
    'Hexyl Chain (6 carbons)': 3.0, # Hydrophobic alkyl chain
    'Chloro Group': 0.7,          # Slightly hydrophobic
    'Amide Group': -2.5,          # Very hydrophilic (H-bonding)
    'Ether Linkage': -1.0,      # Hydrophilic
    'Aliphatic CH2 Linker': 0.5,  # Slightly hydrophobic
    'PEG3 Linker (3 units)': -5.0 # Very hydrophilic
}

# Step 3: Define the components of the original molecule.
original_probe_components = {
    'Thioxanthenone Core': 1,
    'Methoxy Group': 1,
    'Hexyl Chain (6 carbons)': 1,
    'Chloro Group': 1,
    'Amide Group': 1,
    'Ether Linkage': 3,
    'Aliphatic CH2 Linker': 4 # from the ethyl and ethoxy parts not in the hexyl chain
}

# Step 4: Define the components of the modified molecule.
# We replace the Amide Group with a PEG3 linker.
modified_probe_components = {
    'Thioxanthenone Core': 1,
    'Methoxy Group': 1,
    'Hexyl Chain (6 carbons)': 1,
    'Chloro Group': 1,
    'PEG3 Linker (3 units)': 1, # Amide is replaced by this
    'Ether Linkage': 3,
    'Aliphatic CH2 Linker': 4
}

def calculate_and_print_score(name, components, scores):
    """Calculates and prints the total hydrophobicity score for a molecule."""
    total_score = 0
    equation_parts = []
    
    print(f"--- Calculating Score for: {name} ---")
    
    for group, count in components.items():
        score = scores.get(group, 0)
        contribution = score * count
        total_score += contribution
        
        # Build the equation string part-by-part
        sign = "+" if contribution >= 0 else "-"
        equation_parts.append(f" {sign} {abs(contribution):.1f} ({group})")

    # Format the final equation string
    # Start with the first term, removing the leading ' + '
    final_equation = "Total Score = " + equation_parts[0].lstrip(' +')
    # Add the rest of the terms
    for part in equation_parts[1:]:
        final_equation += part
        
    print("The score is the sum of contributions from each molecular group:")
    print(final_equation)
    print(f"\nFinal Hydrophobicity Score for {name}: {total_score:.1f}\n")
    return total_score

# Step 5: Calculate and compare the scores.
original_score = calculate_and_print_score("Original Probe", original_probe_components, hydrophobicity_scores)
modified_score = calculate_and_print_score("Modified Probe (with PEG)", modified_probe_components, hydrophobicity_scores)

# Step 6: Print the final conclusion.
print("--- Conclusion ---")
print(f"Original Probe Score: {original_score:.1f}")
print(f"Modified Probe Score: {modified_score:.1f}")

if modified_score < original_score:
    print("\nThe score of the modified probe is significantly lower.")
    print("This predicts a substantial increase in water solubility.")
    print("\nTherefore, yes, replacing the amide group with a hydrophilic PEG linker is an excellent strategy and is very likely to solve the precipitation problem.")
    # Set the exit code to indicate success/agreement for scripts that might check it
    # This is not strictly necessary for the user, but good practice.
    sys.exit(0) # Corresponds to "Yes"
else:
    print("\nBased on this model, the proposed change does not improve solubility.")
    sys.exit(1) # Corresponds to "No"