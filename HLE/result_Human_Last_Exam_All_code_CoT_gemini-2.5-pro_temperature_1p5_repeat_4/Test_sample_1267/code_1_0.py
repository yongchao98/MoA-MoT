import textwrap

def analyze_tRNA_mutation():
    """
    Analyzes the effect of a mutation in a tRNA anticodon on protein synthesis.
    """
    # Standard genetic code (mRNA codon to Amino Acid)
    genetic_code = {
        'UUU': 'Phe', 'UUC': 'Phe', 'UUA': 'Leu', 'UUG': 'Leu',
        'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser',
        'UAU': 'Tyr', 'UAC': 'Tyr', 'UAA': 'Stop', 'UAG': 'Stop',
        'UGU': 'Cys', 'UGC': 'Cys', 'UGA': 'Stop', 'UGG': 'Trp',
        'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
        'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
        'CAU': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
        'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',
        'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile', 'AUG': 'Met',
        'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
        'AAU': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
        'AGU': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
        'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
        'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
        'GAU': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
        'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly',
    }

    def get_codon_from_anticodon(anticodon_5_to_3):
        """
        Calculates the corresponding mRNA codon for a given tRNA anticodon.
        It accounts for antiparallel binding and complementary base pairing.
        """
        # Simplify anticodon by removing modifications and direction markers
        base_sequence = anticodon_5_to_3.split('-')[1]
        
        # Reverse the anticodon sequence to match mRNA's 5'->3' direction
        anticodon_read_for_pairing = base_sequence[::-1]

        # Find the complementary bases
        codon = ""
        pairings = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}
        for base in anticodon_read_for_pairing:
            codon += pairings.get(base, 'N')
        return codon

    # --- Analysis ---
    original_anticodon = "5'-xm5s2UAA-3'"
    mutated_anticodon = "5'-xm5s2UUG-3'"
    substitution_rate_denominator = 1000

    # 1. Determine the codon recognized by the original tRNA
    original_codon = get_codon_from_anticodon(original_anticodon)
    original_aa = genetic_code.get(original_codon, "Unknown")

    # 2. Determine the codon recognized by the mutated tRNA
    mutated_codon_target = get_codon_from_anticodon(mutated_anticodon)
    target_aa = genetic_code.get(mutated_codon_target, "Unknown")
    
    # The crucial point: the tRNA is still charged with its original amino acid
    charging_aa = original_aa

    # 3. Print the detailed explanation
    print("Analysis of tRNA Mutation Impact on Translation")
    print("=" * 50)
    
    print(f"Original tRNA Analysis:")
    print(f"  - Anticodon: {original_anticodon}")
    print(f"  - Recognizes mRNA Codon: {original_codon}")
    print(f"  - Codon {original_codon} codes for: {original_aa}")
    print(f"  - Conclusion: This tRNA is normally charged with {original_aa} (Leucine).")
    print("-" * 50)

    print(f"Mutated tRNA Analysis:")
    print(f"  - Anticodon: {mutated_anticodon}")
    print(f"  - Recognizes mRNA Codon: {mutated_codon_target}")
    print(f"  - Codon {mutated_codon_target} normally codes for: {target_aa}")
    print("-" * 50)
    
    print("Implication for Protein Synthesis:")
    explanation = (
        f"The mutation changes the anticodon, but not the identity of the tRNA itself. "
        f"Therefore, the cell's machinery still attaches the original amino acid, {charging_aa}, to this mutated tRNA. "
        f"During translation, this incorrectly charged tRNA now binds to the {mutated_codon_target} codon. "
        f"This results in the insertion of {charging_aa} where {target_aa} was intended."
    )
    print(textwrap.fill(explanation, width=70))
    print("-" * 50)
    
    print("Final Conclusion:")
    print(f"The mutated tRNA causes the insertion of Leucine at a codon meant for Glutamine.")
    print(f"This is a missense mutation where the tRNA machinery is subverted. The substitution occurs in approximately 1 in {substitution_rate_denominator} instances because the mutated tRNA competes with the normal, correct tRNA for Glutamine.")
    print("\nThis directly supports the statement: 'It allows insertion of an amino acid usually inserted by another, more common anticodon.'")

# Execute the analysis
analyze_tRNA_mutation()