def calculate_hybrid_stability(rna_sequence):
    """
    Calculates a simplified binding energy for an RNA-DNA hybrid.
    - A-U pairs contribute 2 units of energy.
    - G-C pairs contribute 3 units of energy.
    A higher total energy signifies a more stable hybrid that is harder to dissociate.
    """
    stability_score = 0
    # Create the corresponding DNA template strand (A<->U, G<->C)
    dna_template = ""
    for base in rna_sequence:
        if base == 'U':
            dna_template += 'A'
        elif base == 'A':
            dna_template += 'T'
        elif base == 'G':
            dna_template += 'C'
        elif base == 'C':
            dna_template += 'G'

    # Calculate stability score
    for base in rna_sequence:
        if base == 'U' or base == 'A':
            stability_score += 2
        elif base == 'G' or base == 'C':
            stability_score += 3
            
    return stability_score, dna_template

# The sequence following the 3-4 terminator hairpin
original_attenuator_rna = "UUUUUUUU"
# A hypothetical mutation making the sequence G-C rich
mutated_attenuator_rna = "GCGCGCGC"

# Calculate stability for the original sequence
original_stability, original_dna = calculate_hybrid_stability(original_attenuator_rna)

# Calculate stability for the mutated sequence
mutated_stability, mutated_dna = calculate_hybrid_stability(mutated_attenuator_rna)

print("Analysis of the Rho-Independent Terminator Sequence:")
print("-" * 50)
print(f"The 3-4 terminator hairpin forms, causing RNA polymerase to pause.")
print(f"Termination depends on the stability of the RNA-DNA hybrid that follows.")
print("\nCase 1: Original U-Rich Attenuator")
print(f"RNA Sequence:          {original_attenuator_rna}")
print(f"DNA Template Strand:   {original_dna}")
print(f"Simplified stability score (Sum of 2s for U-A pairs): {' + '.join(['2'] * len(original_attenuator_rna))} = {original_stability}")
print("Result: This low stability allows the RNA to easily dissociate from the DNA, causing transcription to terminate.")

print("\nCase 2: Mutated G-C Rich Sequence")
print(f"RNA Sequence:          {mutated_attenuator_rna}")
print(f"DNA Template Strand:   {mutated_dna}")
print(f"Simplified stability score (Sum of 3s for G-C pairs): {' + '.join(['3'] * len(mutated_attenuator_rna))} = {mutated_stability}")
print("Result: This high stability creates a strong RNA-DNA hybrid that does not dissociate, preventing termination and allowing transcription to continue.")
print("-" * 50)
