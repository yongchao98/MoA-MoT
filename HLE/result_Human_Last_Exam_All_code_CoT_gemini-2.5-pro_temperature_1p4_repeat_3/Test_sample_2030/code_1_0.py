def reverse_complement(dna_sequence):
    """Computes the reverse complement of a DNA sequence."""
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    # Create the complementary strand
    complement_seq = "".join(complement_map.get(base, 'N') for base in dna_sequence)
    # Reverse the complementary strand to get the reverse complement
    reverse_comp_seq = complement_seq[::-1]
    return reverse_comp_seq

def main():
    """
    Solves the DNA oligo design problem based on the provided specifications.
    """
    # Step 1 & 2: Analysis of the 6 reading frames shows Frame 3 has two unique amino acids.
    # Original Sequence: 5' CTT CCC CGC ACA AGT GGT 3'
    # Frame 3 sequence (starting from the 3rd nucleotide 'T'): 5'-TTC CCG CGA CAA GTG-3'
    # Frame 3 translation: Phe(F), Pro(P), Arg(R), Gln(Q), Val(V).
    # The unique amino acids in this frame are Phenylalanine (F) and Glutamine (Q).

    # Step 3: Identify codons and properties.
    # Codon for Phenylalanine (F) is 'TTC'. Phenylalanine is a non-polar amino acid.
    # Codon for Glutamine (Q) is 'CAA'. Glutamine is a polar amino acid.

    # Apply the SNP rules:
    # Rule 1: "changed a polar amino acid to a stop codon"
    # The polar amino acid is Glutamine (Gln, 'CAA'). A single base change to 'TAA' creates a stop codon.
    # Original codon: CAA -> Modified codon: TAA

    # Rule 2: "changed a non-polar amino acid to a cysteine amino acid"
    # The non-polar amino acid is Phenylalanine (Phe, 'TTC'). A single base change to 'TGC' creates a Cysteine codon.
    # Original codon: TTC -> Modified codon: TGC

    # Step 4: Construct the part of the modified sequence that is translated into amino acids.
    # The original codons in the reading frame were: TTC CCG CGA CAA
    # The modified codons become: TGC CCG CGA TAA
    # The translation stops at the 'TAA' codon. Therefore, the sequence that is translated into amino acids is 'TGC CCG CGA'.
    target_sequence = "TGCCCGCGA"
    
    # Step 5: Design the oligo by finding the reverse complement of the target sequence.
    oligo_sequence = reverse_complement(target_sequence)

    # Print the final result in the specified format.
    print(f"Original Frame 3 Codons: TTC CCG CGA CAA")
    print(f"Amino Acid Sequence:      Phe Pro Arg Gln")
    print(f"Modified Frame 3 Codons:  TGC CCG CGA TAA")
    print(f"Modified Amino Acids:     Cys Pro Arg STOP")
    print(f"\nThe target DNA sequence for the oligo (translated part of the modified sequence) is 5'-{target_sequence}-3'.")
    print(f"The DNA sequence of the oligo that binds to this target is its reverse complement.")
    print(f"\nFinal Oligo Sequence: 5' {oligo_sequence} 3'")

if __name__ == "__main__":
    main()
<<<5' TCGCGGGCA 3'>>>