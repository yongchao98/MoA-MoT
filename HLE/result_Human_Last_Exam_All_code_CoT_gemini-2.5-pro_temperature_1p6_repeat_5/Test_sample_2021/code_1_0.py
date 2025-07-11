def solve_rye_problem():
    # Part 1: Genetic Calculation for Drought Resistance

    # Pollination rates
    self_pollination_rate = 0.05
    cross_pollination_rate = 0.95

    # Genotype proportions from self-pollination (WI x WI)
    # Offspring: 1/4 WW, 1/2 WI, 1/4 II
    self_ww = 0.25
    self_wi = 0.50
    self_ii = 0.25

    # Genotype proportions from cross-pollination (WI x WW)
    # Offspring: 1/2 WW, 1/2 WI
    cross_ww = 0.50
    cross_wi = 0.50
    cross_ii = 0.00

    # Total genotype proportions in the next generation
    total_ww = (self_ww * self_pollination_rate) + (cross_ww * cross_pollination_rate)
    total_wi = (self_wi * self_pollination_rate) + (cross_wi * cross_pollination_rate)
    total_ii = (self_ii * self_pollination_rate) + (cross_ii * cross_pollination_rate)

    # The parent plant (WI) is resistant. This implies that having at least one
    # inactive 'I' allele confers resistance. Therefore, WI and II genotypes are resistant.
    resistant_proportion = total_wi + total_ii
    resistant_percentage = resistant_proportion * 100

    # Part 2: Protein Mass Calculation

    # DNA sequence of the 105 nucleotide insertion
    dna_sequence = "gatgaaaaagtgtttaccaaagaactggatcagaccattgaacagctgaacgaatgcggccagctgagcgaaagccaggtgaaaagcctgatgtgcgaagcgacc".upper()

    # Genetic code for translation
    genetic_code = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TGC':'C', 'TGT':'C', 'TGG':'W'
    }

    # Translate DNA to protein
    protein_sequence = ""
    for i in range(0, len(dna_sequence), 3):
        codon = dna_sequence[i:i+3]
        protein_sequence += genetic_code[codon]

    # Monoisotopic residue masses of amino acids (in Daltons)
    aa_masses = {
        'A': 71.03711,  'C': 103.00919, 'D': 115.02694, 'E': 129.04259,
        'F': 147.06841, 'G': 57.02146,  'H': 137.05891, 'I': 113.08406,
        'K': 128.09496, 'L': 113.08406, 'M': 131.04049, 'N': 114.04293,
        'P': 97.05276,  'Q': 128.05858, 'R': 156.10111, 'S': 87.03203,
        'T': 101.04768, 'V': 99.06841,  'W': 186.07931, 'Y': 163.06333
    }
    
    insertion_mass_da = sum(aa_masses[aa] for aa in protein_sequence)
    insertion_mass_kda = insertion_mass_da / 1000

    # Part 3: Synthesize and Print Results
    
    print("--- Analysis Results ---")
    
    print(f"\n1. Theoretical Percentage of Resistant Offspring:")
    print(f"The combined proportion of resistant genotypes (WI and II) is {total_wi:.4f} + {total_ii:.4f} = {resistant_proportion:.4f}.")
    print(f"Theoretically, {resistant_percentage:.2f}% of the offspring should be drought-resistant.")
    
    print(f"\n2. Functional Analysis from Experiment Data:")
    print("The co-expression experiment shows Par22 levels decrease with E3ub-wt (200 units) and increase with E3ub-insert105 (3000 units), compared to control (700 units).")
    print("This indicates that only E3ub-wt is an active ubiquitin ligase that degrades Par22.")
    print("The native mass spectrometry shows a 110kDa complex for 'Par22 plus E3ub-wt', but separate peaks for 'Par22 plus E3ub-insert105', including one at 50kDa (free Par22).")
    print("This indicates that Par22 cannot interact with E3ub-insert105.")
    
    print(f"\n3. Mass Increase from Insertion:")
    print(f"The 105 nucleotide insertion translates to a {len(protein_sequence)}-amino acid peptide.")
    print(f"The calculated mass of this peptide is {insertion_mass_da:.2f} Daltons.")
    print(f"This is an increase in mass for the E3ub protein of approximately {insertion_mass_kda:.1f} kDa.")
    
    print("\n--- Final Conclusion ---")
    print(f"Based on the analysis, the correct statement is: Theoretically, {resistant_percentage:.2f}% of the offspring should be drought-resistant. Only E3ub-wt is an active ubiquitin ligase. Par22 cannot interact with E3ub-insert105. The insertion increases the mass of E3ub by approximately {insertion_mass_kda:.1f} kDA.")


solve_rye_problem()