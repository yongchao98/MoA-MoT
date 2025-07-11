import textwrap

def solve_biology_problem():
    """
    Analyzes the provided biological data to determine the correct statement.
    This involves:
    1. Calculating the percentage of resistant offspring.
    2. Calculating the mass increase from the gene insertion.
    3. Interpreting experimental data on protein activity and interaction.
    """

    # --- Step 1: Genetic Cross Analysis ---
    print("### Step 1: Calculating the Percentage of Resistant Offspring ###")
    print("The parent plant is heterozygous for the resistance allele (Genotype: WI).")
    print("Resistance is dominant, so genotypes WI and II are resistant.")
    print("The surrounding population is wild-type (Genotype: WW).")
    print("\nThere are two pollination scenarios:")

    # Scenario A: Self-pollination (5% probability)
    self_poll_rate = 0.05
    # Offspring from WI x WI cross: 1/4 WW (non-resistant), 1/2 WI (resistant), 1/4 II (resistant)
    # Total resistant fraction from selfing = 1/2 + 1/4 = 3/4 = 0.75
    res_from_selfing = 0.75
    print(f"1. Self-pollination (occurs {self_poll_rate*100}% of the time):")
    print(f"   - Cross: WI x WI -> Offspring are 1/4 WW, 1/2 WI, 1/4 II.")
    print(f"   - Resistant fraction = P(WI) + P(II) = 0.5 + 0.25 = {res_from_selfing}")

    # Scenario B: Cross-pollination (95% probability)
    cross_poll_rate = 0.95
    # Offspring from WI x WW cross: 1/2 WW (non-resistant), 1/2 WI (resistant)
    # Total resistant fraction from crossing = 1/2 = 0.50
    res_from_crossing = 0.50
    print(f"\n2. Cross-pollination with wild-type WW (occurs {cross_poll_rate*100}% of the time):")
    print(f"   - Cross: WI x WW -> Offspring are 1/2 WW, 1/2 WI.")
    print(f"   - Resistant fraction = P(WI) = {res_from_crossing}")

    # Total theoretical percentage of resistant offspring
    total_resistant_fraction = (self_poll_rate * res_from_selfing) + (cross_poll_rate * res_from_crossing)
    total_resistant_percentage = total_resistant_fraction * 100

    print("\nCombining these scenarios to find the total percentage of resistant offspring:")
    print("Final Equation:")
    print(f"({self_poll_rate} * {res_from_selfing}) + ({cross_poll_rate} * {res_from_crossing}) = {total_resistant_fraction}")
    print(f"Total theoretical percentage of resistant offspring = {total_resistant_fraction:.4f} * 100 = {total_resistant_percentage:.2f}%")
    print("-" * 20)

    # --- Step 2: Protein Interaction and Activity Analysis ---
    print("\n### Step 2: Analysis of Protein Function ###")
    print("1. E3 Ligase Activity (Co-expression data):")
    print("   - Par22 + E3ub-wt: Par22 level drops (700 -> 200 units).")
    print("   - Conclusion: E3ub-wt is an active E3 ligase that targets Par22 for degradation.")
    print("   - Par22 + E3ub-insert105: Par22 level increases (700 -> 3000 units), showing no degradation.")
    print("   - Conclusion: E3ub-insert105 is NOT an active E3 ligase towards Par22.")
    print("\n2. Protein-Protein Interaction (Mass Spectrometry data):")
    print("   - Par22 (50kDa) + E3ub-wt (60kDa): A complex of 110kDa is formed (50 + 60).")
    print("   - Conclusion: Par22 interacts with E3ub-wt.")
    print("   - Par22 (50kDa) + E3ub-insert105: Free Par22 at 50kDa is observed.")
    print("   - Conclusion: Par22 does not interact with E3ub-insert105.")
    print("-" * 20)

    # --- Step 3: Mass Calculation of the Insertion ---
    print("\n### Step 3: Calculating the Mass Increase of E3ub Protein ###")
    insertion_dna = "gatgaaaaagtgtttaccaaagaactggatcagaccattgaacagctgaacgaatgcggccagctgagcgaaagccaggtgaaaagcctgatgtgcgaagcgacc"
    
    genetic_code = {
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A', 'UGC': 'C', 'UGU': 'C',
        'GAC': 'D', 'GAU': 'D', 'GAA': 'E', 'GAG': 'E', 'UUC': 'F', 'UUU': 'F',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G', 'CAC': 'H', 'CAU': 'H',
        'AUA': 'I', 'AUC': 'I', 'AUU': 'I', 'AAA': 'K', 'AAG': 'K', 'UUA': 'L',
        'UUG': 'L', 'CUA': 'L', 'CUC': 'L', 'CUG': 'L', 'CUU': 'L', 'AUG': 'M',
        'AAC': 'N', 'AAU': 'N', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P',
        'CAA': 'Q', 'CAG': 'Q', 'AGA': 'R', 'AGG': 'R', 'CGA': 'R', 'CGC': 'R',
        'CGG': 'R', 'CGU': 'R', 'AGC': 'S', 'AGU': 'S', 'UCA': 'S', 'UCC': 'S',
        'UCG': 'S', 'UCU': 'S', 'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACU': 'T',
        'GUA': 'V', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V', 'UGG': 'W', 'UAC': 'Y',
        'UAU': 'Y', 'UAA': '*', 'UAG': '*', 'UGA': '*'
    }
    
    # Average molecular weight of amino acids in Daltons (Da)
    aa_weights = {
        'A': 89.09, 'C': 121.16, 'D': 133.10, 'E': 147.13, 'F': 165.19,
        'G': 75.07, 'H': 155.16, 'I': 131.17, 'K': 146.19, 'L': 131.17,
        'M': 149.21, 'N': 132.12, 'P': 115.13, 'Q': 146.15, 'R': 174.20,
        'S': 105.09, 'T': 119.12, 'V': 117.15, 'W': 204.23, 'Y': 181.19
    }
    
    mass_of_water = 18.015
    
    # Convert DNA to RNA by replacing 't' with 'u'
    rna_sequence = insertion_dna.upper().replace('T', 'U')
    
    codons = textwrap.wrap(rna_sequence, 3)
    protein_sequence = "".join([genetic_code.get(codon, '?') for codon in codons])
    
    total_mass = sum(aa_weights.get(aa, 0) for aa in protein_sequence)
    # Subtract the mass of water for each peptide bond
    peptide_bonds = len(protein_sequence) - 1
    final_mass_da = total_mass - (peptide_bonds * mass_of_water)
    final_mass_kda = final_mass_da / 1000
    
    print(f"The 105 nucleotide insertion codes for {len(protein_sequence)} amino acids.")
    print(f"The calculated mass of this peptide is {final_mass_da:.2f} Da.")
    print(f"This corresponds to an increase of approximately {final_mass_kda:.1f} kDa.")
    print("-" * 20)

    # --- Step 4: Final Conclusion ---
    print("\n### Final Conclusion ###")
    print("Based on the analysis:")
    print(f"1. The theoretical percentage of drought-resistant offspring is {total_resistant_percentage:.2f}%.")
    print("2. Only the E3ub-wt variant is an active ubiquitin ligase.")
    print("3. Par22 cannot interact with E3ub-insert105 (only with E3ub-wt).")
    print(f"4. The insertion increases the mass of the E3ub protein by approximately {final_mass_kda:.1f} kDa (~4.0 kDa).")
    print("\nThis matches all points in answer choice J.")


solve_biology_problem()
<<<J>>>