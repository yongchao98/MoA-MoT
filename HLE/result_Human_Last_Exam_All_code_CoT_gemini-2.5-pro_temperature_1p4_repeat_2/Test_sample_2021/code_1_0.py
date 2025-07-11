def solve_biology_problem():
    """
    Performs calculations for offspring probability and protein mass increase.
    """

    # --- Part 1: Calculate the probability of resistant offspring ---

    # Pollination rates
    self_pollination_rate = 0.05
    cross_pollination_rate = 0.95

    # Genotype probabilities from self-pollination (WI x WI)
    # Resistant genotypes are WI (heterozygous) and II (homozygous for insertion)
    prob_wi_from_self = 0.50
    prob_ii_from_self = 0.25
    prob_resistant_from_self = prob_wi_from_self + prob_ii_from_self

    # Genotype probabilities from cross-pollination (WI x WW)
    # The only resistant genotype is WI
    prob_wi_from_cross = 0.50
    prob_resistant_from_cross = prob_wi_from_cross

    # Total probability of resistant offspring
    total_prob_resistant = (prob_resistant_from_self * self_pollination_rate) + \
                           (prob_resistant_from_cross * cross_pollination_rate)

    print("--- Offspring Resistance Calculation ---")
    print(f"The calculation for the percentage of resistant offspring is:")
    print(f"(P(WI|self) + P(II|self)) * P(self) + P(WI|cross) * P(cross)")
    print(f"= ({prob_wi_from_self} + {prob_ii_from_self}) * {self_pollination_rate} + {prob_wi_from_cross} * {cross_pollination_rate}")
    print(f"= {prob_resistant_from_self} * {self_pollination_rate} + {prob_resistant_from_cross} * {cross_pollination_rate}")
    print(f"= {prob_resistant_from_self * self_pollination_rate} + {prob_resistant_from_cross * cross_pollination_rate}")
    print(f"= {total_prob_resistant}")
    print(f"Theoretically, {total_prob_resistant:.2%} of the offspring should be drought-resistant.\n")

    # --- Part 2: Calculate the mass increase from the DNA insertion ---

    dna_insertion_seq = "gatgaaaaagtgtttaccaaagaactggatcagaccattgaacagctgaacgaatgcggccagctgagcgaaagccaggtgaaaagcctgatgtgcgaagcgacc"
    
    codon_table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }

    # Average mass of amino acid residues in Daltons (Da)
    residue_mass = {
        'A': 71.0788, 'R': 156.1875, 'N': 114.1038, 'D': 115.0886, 'C': 103.1388,
        'E': 129.1155, 'Q': 128.1307, 'G': 57.0519, 'H': 137.1411, 'I': 113.1594,
        'L': 113.1594, 'K': 128.1741, 'M': 131.1926, 'F': 147.1766, 'P': 97.1167,
        'S': 87.0782, 'T': 101.1051, 'W': 186.2132, 'Y': 163.1760, 'V': 99.1326
    }

    num_nucleotides = len(dna_insertion_seq)
    num_amino_acids = num_nucleotides // 3
    
    protein_sequence = ""
    for i in range(0, num_nucleotides, 3):
        codon = dna_insertion_seq[i:i+3]
        protein_sequence += codon_table[codon]
    
    total_mass_increase_da = sum(residue_mass[aa] for aa in protein_sequence)
    total_mass_increase_kda = total_mass_increase_da / 1000

    print("--- Mass Increase Calculation ---")
    print(f"The insertion sequence is {num_nucleotides} nucleotides long.")
    print(f"This translates to {num_nucleotides} / 3 = {num_amino_acids} amino acids.")
    print(f"The calculated mass increase is {total_mass_increase_da:.2f} Da, which is approximately {total_mass_increase_kda:.1f} kDa.\n")

    # --- Part 3: Final Conclusion ---
    print("--- Final Conclusion ---")
    print("Based on the data:")
    print(f"Theoretically, {total_prob_resistant:.2%} of the offspring should be drought-resistant.")
    print("Only E3ub-wt is an active ubiquitin ligase.")
    print("Par22 cannot interact with E3ub-insert105.")
    print(f"The insertion increases the mass of E3ub by approximately {total_mass_increase_kda:.1f} kDA.")

solve_biology_problem()
<<<J>>>