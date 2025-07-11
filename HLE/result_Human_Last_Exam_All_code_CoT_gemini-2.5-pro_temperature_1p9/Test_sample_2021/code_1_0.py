import sys
# Redirect print to stdout to ensure output is displayed
# This is a good practice for environments where stdout might be captured.
sys.stdout = sys.__stdout__

def solve_rye_genetics():
    """
    Analyzes the provided biological data to find the correct conclusion.
    """
    
    # Part 1: Logic Derivations from the problem description
    # 1a. Protein Function from co-expression data:
    # - E3ub-wt reduces Par22 levels (700 -> 200), so it is an ACTIVE E3 ubiquitin ligase.
    # - E3ub-insert105 increases Par22 levels (700 -> 3000), so it is INACTIVE.
    
    # 1b. Protein Interaction from mass spectrometry data:
    # - E3ub-wt (60kDa) + Par22 (50kDa) form a 110kDa complex, so they INTERACT.
    # - E3ub-insert105 + Par22 show separate peaks, so they DO NOT interact.
    
    # 1c. Resistance Mechanism:
    # - The parent plant is heterozygous (E3ub-wt/E3ub-insert105) and resistant.
    # - Resistance is linked to the INACTIVE insert allele, which prevents Par22 degradation.
    # - Therefore, genotypes with at least one insert allele (heterozygous or homozygous for the insertion) are resistant.
    
    # Part 2: Calculation of the insertion's mass contribution
    insertion_dna = "gatgaaaaagtgtttaccaaagaactggatcagaccattgaacagctgaacgaatgcggccagctgagcgaaagccaggtgaaaagcctgatgtgcgaagcgacc"
    nucleotide_count = len(insertion_dna)
    amino_acid_count = nucleotide_count / 3
    
    # Using average residue mass to approximate. For a more precise calculation:
    # DNA to Protein sequence: D E K V F T K E L D Q T I E Q L N E C G Q L S E S Q V K S L M C E A T
    residue_mass = {
        'A': 71.08, 'C': 103.14, 'D': 115.09, 'E': 129.12, 'F': 147.18,
        'G': 57.05, 'I': 113.16, 'K': 128.17, 'L': 113.16, 'M': 131.20,
        'N': 114.10, 'Q': 128.13, 'S': 87.08, 'T': 101.11, 'V': 99.13,
        'R': 156.19, 'H': 137.14, 'P': 97.12, 'W': 186.21, 'Y': 163.18
    }
    protein_sequence = "DEKVFTKELDQTIEQLNECGQLSESQVKSLMCEAT"
    mass_increase_daltons = sum(residue_mass.get(aa, 0) for aa in protein_sequence)
    mass_increase_kda = mass_increase_daltons / 1000

    # Part 3: Calculation of the resistant offspring percentage
    self_pollination_rate = 0.05
    cross_pollination_rate = 0.95
    
    # Parent plant is heterozygous (let's use Ee notation).
    # Resistance is conferred by the 'e' allele (insertion). Genotypes Ee and ee are resistant.
    
    # Offspring from self-pollination (Ee x Ee): 1/4 EE, 2/4 Ee, 1/4 ee.
    # Resistant fraction = Ee + ee = 2/4 + 1/4 = 3/4
    resistant_fraction_selfing = 0.75
    
    # Offspring from cross-pollination with general population (Ee x EE): 1/2 EE, 1/2 Ee
    # Resistant fraction = Ee = 1/2
    resistant_fraction_crossing = 0.50
    
    # Total probability of resistant offspring
    total_resistant_fraction = (self_pollination_rate * resistant_fraction_selfing) + (cross_pollination_rate * resistant_fraction_crossing)
    total_resistant_percentage = total_resistant_fraction * 100
    
    # Print summary and results
    print("Summary of Key Findings:")
    print("1. Activity: Only E3ub-wt is an active ubiquitin ligase.")
    print("2. Interaction: Par22 cannot interact with E3ub-insert105.")
    print(f"3. Mass Increase: The insertion adds ~{mass_increase_kda:.2f} kDa to the protein mass.")
    
    print("\nCalculation for Percentage of Resistant Offspring:")
    print(f"Parent genotype is heterozygous and resistant.")
    print(f"Resistant genotypes in offspring are heterozygous and homozygous for the insertion.")
    
    print("\nThe final equation is based on pollination rates:")
    print(f"Total Resistant % = (Self-Pollination Rate * Resistant Fraction from Selfing) + (Cross-Pollination Rate * Resistant Fraction from Crossing)")
    print(f"Total Resistant Fraction = ({self_pollination_rate} * {resistant_fraction_selfing}) + ({cross_pollination_rate} * {resistant_fraction_crossing})")
    print(f"Total Resistant Fraction = {total_resistant_fraction}")
    print(f"Result: The predicted percentage of drought-resistant offspring is {total_resistant_percentage:.2f}%.")
    
    print("\nConclusion: The statement that aligns with all findings is J.")

# Execute the analysis
solve_rye_genetics()