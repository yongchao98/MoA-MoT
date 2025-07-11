import sys

def solve_rye_genetics_problem():
    """
    Analyzes the provided biological data to determine protein function,
    interaction, mass changes, and calculates the theoretical percentage of
    drought-resistant offspring.
    """
    
    # --- Data from the problem description ---
    
    # Part 1 & 2: Functional and interaction data analysis
    # This part is based on logical deduction from the text.
    print("Step 1 & 2: Analysis of Protein Function and Interaction")
    print("----------------------------------------------------------")
    print("Functional Analysis (Co-expression experiment):")
    print("- Control (Par22 alone): 700 units.")
    print("- Par22 + E3ub-wt: 200 units. The level of Par22 is reduced, so E3ub-wt degrades Par22.")
    print("  Conclusion: E3ub-wt is an ACTIVE E3 ubiquitin ligase.")
    print("- Par22 + E3ub-insert105: 3000 units. The level of Par22 is increased, so E3ub-insert105 does not degrade Par22.")
    print("  Conclusion: E3ub-insert105 is an INACTIVE E3 ubiquitin ligase.")
    print("\nInteraction Analysis (Native Mass Spectrometry):")
    print("- Par22 (50 kDa) + E3ub-wt (60 kDa) -> 110 kDa complex.")
    print("  Conclusion: E3ub-wt INTERACTS with Par22.")
    print("- Par22 (50 kDa) + E3ub-insert105 -> Peak at 50 kDa (free Par22).")
    print("  Conclusion: E3ub-insert105 does NOT interact with Par22.")
    
    # Part 3: Mass calculation
    print("\nStep 3: Calculation of Mass Increase from Insertion")
    print("----------------------------------------------------")
    dna_sequence = "gatgaaaaagtgtttaccaaagaactggatcagaccattgaacagctgaacgaatgcggccagctgagcgaaagccaggtgaaaagcctgatgtgcgaagcgacc"
    
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
    
    aa_residue_mass = {
        'A': 71.08, 'R': 156.19, 'N': 114.11, 'D': 115.09, 'C': 103.14, 'E': 129.12, 
        'Q': 128.14, 'G': 57.05, 'H': 137.14, 'I': 113.16, 'L': 113.16, 'K': 128.17, 
        'M': 131.20, 'F': 147.18, 'P': 97.12, 'S': 87.08, 'T': 101.11, 'W': 186.21, 
        'Y': 163.18, 'V': 99.13
    }
    H2O_mass = 18.015

    if len(dna_sequence) % 3 != 0:
        print("Error: DNA sequence length is not a multiple of 3.", file=sys.stderr)
        return
        
    protein_sequence = ""
    for i in range(0, len(dna_sequence), 3):
        codon = dna_sequence[i:i+3].upper()
        protein_sequence += codon_table.get(codon, '?')

    peptide_mass = H2O_mass + sum(aa_residue_mass.get(aa, 0) for aa in protein_sequence)
    mass_in_kda = peptide_mass / 1000

    print(f"The insertion DNA has {len(dna_sequence)} nucleotides.")
    print(f"This translates into a peptide of {len(protein_sequence)} amino acids.")
    print(f"The calculated mass of this inserted peptide is {peptide_mass:.2f} Daltons.")
    print(f"This is an increase of {mass_in_kda:.2f} kDa, which is approximately 4.0 kDa.")
    
    # Part 4: Offspring calculation
    print("\nStep 4: Calculation of Resistant Offspring Percentage")
    print("-------------------------------------------------------")
    self_pollination_rate = 0.05
    cross_pollination_rate = 0.95
    
    # Parent is heterozygous (Rr), resistance is dominant (R = insert105 allele).
    # Self-pollination (Rr x Rr) gives 3/4 resistant offspring (RR, Rr, rR).
    prob_resistant_selfing = 0.75
    
    # Cross-pollination (Rr x rr) gives 1/2 resistant offspring (Rr).
    prob_resistant_crossing = 0.50
    
    total_prob_resistant = (self_pollination_rate * prob_resistant_selfing) + \
                           (cross_pollination_rate * prob_resistant_crossing)
                           
    total_percent_resistant = total_prob_resistant * 100
    
    print("The plant is heterozygous. Resistance is conferred by the inactive E3ub-insert105 allele.")
    print("Contribution from self-pollination (5% of the time):")
    print(f"  Rate: {self_pollination_rate}, Fraction of resistant offspring: {prob_resistant_selfing}")
    print("Contribution from cross-pollination with wild-type (95% of the time):")
    print(f"  Rate: {cross_pollination_rate}, Fraction of resistant offspring: {prob_resistant_crossing}")
    
    print("\nThe final equation for the total proportion of resistant offspring is:")
    print(f"({self_pollination_rate} * {prob_resistant_selfing}) + ({cross_pollination_rate} * {prob_resistant_crossing}) = {total_prob_resistant:.4f}")
    print(f"This means {total_percent_resistant:.2f}% of the offspring will be resistant.")

    print("\nStep 5: Final Conclusion")
    print("--------------------------")
    print("The results match option J: 'Theoretically, 51.25% of the offspring should be drought-resistant. Only E3ub-wt is an active ubiquitin ligase. Par22 cannot interact with E3ub-insert105. The insertion increases the mass of E3ub by approximately 4.0kDA.'")

solve_rye_genetics_problem()
<<<J>>>