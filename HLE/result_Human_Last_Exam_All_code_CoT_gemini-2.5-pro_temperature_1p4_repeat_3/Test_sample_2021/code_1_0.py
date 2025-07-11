import math

def solve_biology_problem():
    """
    This script performs the calculations needed to solve the given problem by:
    1. Calculating the theoretical frequency of drought-resistant offspring.
    2. Calculating the mass increase of the E3ub protein due to the nucleotide insertion.
    3. Printing a summary statement that corresponds to the correct answer choice.
    """

    # --- Part 1: Calculate the theoretical frequency of resistant offspring ---
    self_pollination_rate = 0.05
    cross_pollination_rate = 0.95
    prop_resistant_selfing = 0.75  # Rr x Rr -> 75% resistant (RR, Rr)
    prop_resistant_crossing = 0.50  # Rr x rr -> 50% resistant (Rr)

    # The equation for total resistant proportion:
    total_resistant_prop = (self_pollination_rate * prop_resistant_selfing) + (cross_pollination_rate * prop_resistant_crossing)
    total_resistant_percentage = total_resistant_prop * 100

    print("--- Offspring Resistance Calculation ---")
    print(f"The equation for the frequency of resistant offspring is:")
    print(f"({self_pollination_rate} * {prop_resistant_selfing}) + ({cross_pollination_rate} * {prop_resistant_crossing}) = {total_resistant_prop}")
    print(f"Resulting percentage: {total_resistant_percentage:.2f}%\n")


    # --- Part 2: Calculate the mass increase from the nucleotide insertion ---
    dna_insertion = "gatgaaaaagtgtttaccaaagaactggatcagaccattgaacagctgaacgaatgcggccagctgagcgaaagccaggtgaaaagcctgatgtgcgaagcgacc"
    
    # Standard genetic code and monoisotopic residue masses (in Daltons)
    codon_map = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',}
    aa_mass = {
        'A': 71.03711, 'R': 156.10111, 'N': 114.04293, 'D': 115.02694, 'C': 103.00919, 
        'E': 129.04259, 'Q': 128.05858, 'G': 57.02146, 'H': 137.05891, 'I': 113.08406, 
        'L': 113.08406, 'K': 128.09496, 'M': 131.04049, 'F': 147.06841, 'P': 97.05276, 
        'S': 87.03203, 'T': 101.04768, 'W': 186.07931, 'Y': 163.06333, 'V': 99.06841, '_': 0}
    mass_H2O = 18.01528

    peptide_sequence = ""
    for i in range(0, len(dna_insertion), 3):
        codon = dna_insertion[i:i+3].upper()
        peptide_sequence += codon_map.get(codon, '?')

    sum_of_residue_masses = sum(aa_mass.get(aa, 0) for aa in peptide_sequence)
    peptide_mass_daltons = sum_of_residue_masses + mass_H2O
    peptide_mass_kda = peptide_mass_daltons / 1000

    print("--- Mass Increase Calculation ---")
    print(f"The 105 nucleotide insertion codes for {len(peptide_sequence)} amino acids.")
    print(f"The equation for the peptide's mass is: Sum of residue masses + Mass of H2O")
    print(f"{sum_of_residue_masses:.2f} Da + {mass_H2O:.2f} Da = {peptide_mass_daltons:.2f} Da")
    print(f"Mass increase in kDa: {peptide_mass_kda:.1f} kDA\n") # Rounded to one decimal place for "approximately 4.0"

    # --- Part 3: Final Conclusion ---
    print("--- Final Conclusion ---")
    print(
        f"Theoretically, {total_resistant_percentage:.2f}% of the offspring should be drought-resistant. "
        f"Only E3ub-wt is an active ubiquitin ligase. "
        f"Par22 cannot interact with E3ub-insert105. "
        f"The insertion increases the mass of E3ub by approximately {peptide_mass_kda:.1f}kDA."
    )

solve_biology_problem()