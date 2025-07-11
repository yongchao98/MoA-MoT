def reverse_complement(seq):
    """Computes the reverse complement of an RNA sequence."""
    complement_map = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}
    return "".join(complement_map[base] for base in reversed(seq))

def main():
    """
    Analyzes the effect of a tRNA anticodon mutation.
    """
    genetic_code = {
        'UUA': 'Leucine (Leu)', 'UUG': 'Leucine (Leu)',
        'CUU': 'Leucine (Leu)', 'CUC': 'Leucine (Leu)',
        'CUA': 'Leucine (Leu)', 'CUG': 'Leucine (Leu)',
        'CAA': 'Glutamine (Gln)', 'CAG': 'Glutamine (Gln)',
        # Add other codons if needed for a more complete model
    }

    # --- Original tRNA ---
    original_anticodon = 'UAA' # Core sequence from 5'-xm5s2UAA-3'
    original_tRNA_carries = 'Leucine'
    # The anticodon pairs with the mRNA codon
    original_target_codon = reverse_complement(original_anticodon)
    original_codon_codes_for = genetic_code.get(original_target_codon, 'Unknown')

    print("--- Original State ---")
    print(f"Original tRNA anticodon (5'->3'): {original_anticodon}")
    print(f"This tRNA carries the amino acid: {original_tRNA_carries}")
    print(f"It recognizes the mRNA codon (5'->3'): {original_target_codon}")
    print(f"The codon {original_target_codon} correctly codes for: {original_codon_codes_for}\n")

    # --- Mutated tRNA ---
    mutated_anticodon = 'UUG' # Core sequence from 5'-xm5s2UUG-3'
    # The mutated tRNA still carries Leucine
    mutated_tRNA_carries = 'Leucine'
    # The mutated anticodon recognizes a new mRNA codon
    new_target_codon = reverse_complement(mutated_anticodon)
    new_codon_normally_codes_for = genetic_code.get(new_target_codon, 'Unknown')

    print("--- Mutated State ---")
    print(f"Mutated tRNA anticodon (5'->3'): {mutated_anticodon}")
    print(f"This mutated tRNA still carries: {mutated_tRNA_carries}")
    print(f"It now mistakenly recognizes the mRNA codon (5'->3'): {new_target_codon}")
    print(f"The codon {new_target_codon} should normally code for: {new_codon_normally_codes_for}\n")

    # --- Implication ---
    print("--- Implication ---")
    print(f"The mutation causes {mutated_tRNA_carries} to be inserted where {new_codon_normally_codes_for} should be.")
    print("This happens when the ribosome encounters the codon {} in an mRNA sequence.".format(new_target_codon))
    print("This scenario matches answer choice C: It allows insertion of an amino acid (Leucine) usually inserted by another, more common anticodon (the one for Glutamine).")


if __name__ == "__main__":
    main()

<<<C>>>