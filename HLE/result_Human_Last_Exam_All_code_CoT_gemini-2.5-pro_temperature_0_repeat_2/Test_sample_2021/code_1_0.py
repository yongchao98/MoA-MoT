import textwrap

def solve_biology_problem():
    # Step 1: Calculate the percentage of resistant offspring
    self_pollination_rate = 0.05
    cross_pollination_rate = 0.95

    # Genotype frequencies from self-pollination (wt/ins x wt/ins)
    self_wt_ins = 0.50
    self_ins_ins = 0.25

    # Genotype frequencies from cross-pollination (wt/ins x wt/wt)
    cross_wt_ins = 0.50
    cross_ins_ins = 0.0

    # Total frequency of resistant genotypes
    total_wt_ins = (self_pollination_rate * self_wt_ins) + (cross_pollination_rate * cross_wt_ins)
    total_ins_ins = (self_pollination_rate * self_ins_ins) + (cross_pollination_rate * cross_ins_ins)

    total_resistant_freq = total_wt_ins + total_ins_ins
    total_resistant_percent = total_resistant_freq * 100

    # Step 2: Calculate the mass increase from the insertion
    insertion_dna = "gatgaaaaagtgtttaccaaagaactggatcagaccattgaacagctgaacgaatgcggccagctgagcgaaagccaggtgaaaagcctgatgtgcgaagcgacc"

    genetic_code = {
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R', 'AGA': 'R', 'AGG': 'R',
        'AAC': 'N', 'AAT': 'N',
        'GAC': 'D', 'GAT': 'D',
        'TGC': 'C', 'TGT': 'C',
        'GAA': 'E', 'GAG': 'E',
        'CAA': 'Q', 'CAG': 'Q',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'CAC': 'H', 'CAT': 'H',
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L', 'TTA': 'L', 'TTG': 'L',
        'AAA': 'K', 'AAG': 'K',
        'ATG': 'M',
        'TTC': 'F', 'TTT': 'F',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'AGC': 'S', 'AGT': 'S', 'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'TGG': 'W',
        'TAC': 'Y', 'TAT': 'Y',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'TAA': '_', 'TAG': '_', 'TGA': '_'
    }

    # Average isotopic masses of amino acid residues in Daltons (Da)
    residue_masses = {
        'A': 71.0788, 'R': 156.1875, 'N': 114.1038, 'D': 115.0886, 'C': 103.1388,
        'E': 129.1155, 'Q': 128.1307, 'G': 57.0519, 'H': 137.1411, 'I': 113.1594,
        'L': 113.1594, 'K': 128.1741, 'M': 131.1926, 'F': 147.1766, 'P': 97.1167,
        'S': 87.0782, 'T': 101.1051, 'W': 186.2132, 'Y': 163.1760, 'V': 99.1326
    }

    peptide_sequence = ""
    for i in range(0, len(insertion_dna), 3):
        codon = insertion_dna[i:i+3].upper()
        peptide_sequence += genetic_code[codon]

    mass_increase_daltons = sum(residue_masses[aa] for aa in peptide_sequence)
    mass_increase_kda = mass_increase_daltons / 1000

    # Step 3: Print the reasoning and final conclusion
    print("Step-by-step Analysis:")
    print("1. Calculation of Resistant Offspring Percentage:")
    print(f"   - The parent plant (wt/ins) has a {self_pollination_rate*100}% self-pollination rate and a {cross_pollination_rate*100}% cross-pollination rate with wild-type (wt/wt) plants.")
    print(f"   - Contribution from selfing to resistant genotypes (wt/ins + ins/ins): ({self_pollination_rate} * {self_wt_ins}) + ({self_pollination_rate} * {self_ins_ins}) = {self_pollination_rate * (self_wt_ins + self_ins_ins)}")
    print(f"   - Contribution from crossing to resistant genotypes (wt/ins): {cross_pollination_rate} * {cross_wt_ins} = {cross_pollination_rate * cross_wt_ins}")
    print(f"   - Total resistant offspring = {total_wt_ins:.4f} (from wt/ins) + {total_ins_ins:.4f} (from ins/ins) = {total_resistant_freq:.4f}")
    print(f"   - Therefore, the theoretical percentage of drought-resistant offspring is {total_resistant_percent:.2f}%.")
    print("\n2. Analysis of Protein Function:")
    print("   - E3 Ligase Activity: The co-expression data shows Par22 levels drop with E3ub-wt but rise with E3ub-insert105. This indicates that only E3ub-wt is an active ubiquitin ligase targeting Par22.")
    print("   - Protein Interaction: The native mass spectrometry shows a complex (110kDa) for Par22 + E3ub-wt, but no complex for Par22 + E3ub-insert105. This indicates that Par22 cannot interact with E3ub-insert105.")
    print("\n3. Calculation of Protein Mass Increase:")
    print(f"   - The 105 nucleotide insertion translates to {len(peptide_sequence)} amino acids.")
    print(f"   - The calculated mass of the inserted peptide is {mass_increase_daltons:.2f} Da, which is approximately {mass_increase_kda:.2f} kDa.")
    print("\nConclusion:")
    print("Combining these findings, the correct statement is:")
    print(f"Theoretically, {total_resistant_percent:.2f}% of the offspring should be drought-resistant. Only E3ub-wt is an active ubiquitin ligase. Par22 cannot interact with E3ub-insert105. The insertion increases the mass of E3ub by approximately {mass_increase_kda:.1f}kDA.")

solve_biology_problem()
<<<J>>>