def solve_oligo_design():
    """
    Solves the oligo design puzzle based on the provided constraints.
    """

    # Step 1: Identify the target DNA sequence based on the problem's logic.
    # After analyzing all 6 reading frames, the only one that fits the SNP mutation criteria
    # (Polar -> Stop, Non-polar -> Cys, no initial Cys) is the +3 frame, which translates to
    # Ser-Arg-Ala-Gln-Trp. The SNPs mutate Gln to a Stop codon and Trp to Cys.
    # The oligo should bind to the sequence that is still translated, which are the codons
    # before the new Stop codon.
    # Codons for Ser-Arg-Ala: TCC, CGC, GCA
    target_dna = "TCCCGCGCA"

    # Step 2: Calculate the reverse complement of the target DNA.
    # This will be the sequence of the oligo.
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    
    # First, find the complement strand
    try:
        complement_dna = "".join([complement_map[base] for base in target_dna])
    except KeyError as e:
        print(f"Error: Invalid character {e} in DNA sequence.")
        return

    # Then, reverse the complement strand to get the reverse complement
    reverse_complement = complement_dna[::-1]

    # Step 3: Print the final oligo sequence with its orientation.
    print("The DNA sequence of the oligo that binds to the translated portion of the modified sequence is:")
    print(f"5' {reverse_complement} 3'")

solve_oligo_design()