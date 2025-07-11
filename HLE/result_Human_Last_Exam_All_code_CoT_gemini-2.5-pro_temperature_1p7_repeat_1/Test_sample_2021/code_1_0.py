import sys
# Python 2 and 3 compatibility
if sys.version_info[0] < 3:
    # Redefine print for Python 2
    _print = print
    def print(*args, **kwargs):
        _print(*args)
        sys.stdout.flush()

def solve_biology_problem():
    """
    This function breaks down the problem, performs the necessary calculations,
    and explains the reasoning to arrive at the correct answer.
    """
    
    print("Step 1: Calculate the theoretical percentage of resistant offspring.")
    
    # Pollination rates
    self_pollination_rate = 0.05  # 5%
    cross_pollination_rate = 0.95 # 95%
    
    # Genetic assumptions based on the answer choices matching 51.25%:
    # The parent plant is heterozygous (wt/insert).
    # Resistance is dominant, meaning both heterozygous (wt/insert) and homozygous mutant (insert/insert) are resistant.
    
    # Probability of resistant offspring from self-pollination (wt/insert x wt/insert)
    # Offspring genotypes: 1/4 wt/wt (non-resistant), 1/2 wt/insert (resistant), 1/4 insert/insert (resistant)
    # Total resistant fraction = 1/2 + 1/4 = 3/4
    resistance_from_selfing = 0.75
    
    # Probability of resistant offspring from cross-pollination (wt/insert x wt/wt)
    # The general population is assumed to be wt/wt (non-resistant).
    # Offspring genotypes: 1/2 wt/wt (non-resistant), 1/2 wt/insert (resistant)
    # Total resistant fraction = 1/2
    resistance_from_crossing = 0.50
    
    # Total percentage of resistant offspring
    total_resistant_percentage = (self_pollination_rate * resistance_from_selfing + cross_pollination_rate * resistance_from_crossing) * 100
    
    print("The resistance allele (insert) acts as a dominant allele.")
    print("Offspring from self-pollination (5% of total) have a 3/4 chance of being resistant.")
    print("Offspring from cross-pollination (95% of total) have a 1/2 chance of being resistant.")
    print("Total resistant offspring percentage = ({} * {}) + ({} * {}) = {:.2f}%".format(self_pollination_rate, resistance_from_selfing, cross_pollination_rate, resistance_from_crossing, total_resistant_percentage))
    
    print("\nStep 2: Interpret the biochemical data.")
    print("Densitometry data shows that Par22 protein levels decrease only when co-expressed with E3ub-wt.")
    print("This indicates that only E3ub-wt is an active E3 ubiquitin ligase that targets Par22 for degradation.")
    print("Mass spectrometry data shows a complex (110kDa) forming between Par22 (50kDa) and E3ub-wt (60kDa).")
    print("When Par22 is mixed with E3ub-insert105, free Par22 (50kDa) is detected, showing a lack of stable interaction.")
    print("Conclusion: Par22 cannot interact with E3ub-insert105.")
    
    print("\nStep 3: Calculate the mass increase from the insertion.")
    
    nucleotide_insertion_length = 105 # base pairs
    # An in-frame insertion means the number of nucleotides is a multiple of 3.
    amino_acids_added = nucleotide_insertion_length / 3.0
    
    # Average molecular weight of an amino acid is approximately 110 Daltons (Da).
    avg_amino_acid_mass_da = 110
    
    mass_increase_da = amino_acids_added * avg_amino_acid_mass_da
    mass_increase_kda = mass_increase_da / 1000.0
    
    print("The insertion has {} nucleotides.".format(nucleotide_insertion_length))
    print("This corresponds to {} / 3 = {} amino acids.".format(nucleotide_insertion_length, int(amino_acids_added)))
    print("Mass increase = {} amino acids * ~{} Da/amino acid = {} Da, which is approximately {:.1f} kDa.".format(int(amino_acids_added), avg_amino_acid_mass_da, int(mass_increase_da), mass_increase_kda))
    
    print("\nConclusion: Based on all the evidence, the correct statement must include:")
    print("- Theoretical resistance of 51.25%")
    print("- Only E3ub-wt is an active ubiquitin ligase.")
    print("- Par22 cannot interact with E3ub-insert105.")
    print("- Mass increase is approximately 4.0 kDa.")
    print("This matches option J.")

solve_biology_problem()