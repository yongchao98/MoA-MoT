import math

def solve_biology_problem():
    """
    This function analyzes the provided biological data to determine the correct answer choice.
    It breaks down the analysis into four parts:
    1. E3ub ligase activity
    2. Protein interaction
    3. Mass increase from the insertion
    4. Genetic calculation of offspring resistance
    """
    
    print("Step-by-step analysis of the experimental data:")

    # 1. Analyze E3ub ligase activity
    par22_alone = 700
    par22_with_wt = 200
    par22_with_insert = 3000
    print("\n1. E3ub Ligase Activity Analysis:")
    print(f" - Par22 expressed alone (control): {par22_alone} units.")
    print(f" - Par22 co-expressed with E3ub-wt: {par22_with_wt} units.")
    print("   The significant decrease in Par22 level indicates that E3ub-wt is actively targeting Par22 for degradation.")
    print(f" - Par22 co-expressed with E3ub-insert105: {par22_with_insert} units.")
    print("   The high level of Par22 indicates it is not being degraded by E3ub-insert105.")
    print("Conclusion: Only E3ub-wt is an active ubiquitin ligase towards Par22.")

    # 2. Analyze protein interaction
    par22_mass = 50
    e3ub_wt_mass = 60
    complex_mass = 110
    print("\n2. Protein Interaction Analysis (Native Mass Spectrometry):")
    print(f" - Par22 ({par22_mass} kDa) + E3ub-wt ({e3ub_wt_mass} kDa) result in a peak at {complex_mass} kDa.")
    print(f"   Since {par22_mass} + {e3ub_wt_mass} = {complex_mass}, this indicates a stable 1:1 complex is formed, confirming interaction.")
    print(f" - Par22 ({par22_mass} kDa) + E3ub-insert105 results in a peak at {par22_mass} kDa (free Par22).")
    print("   The presence of free Par22 and the absence of a complex peak indicates that Par22 does not interact with E3ub-insert105.")
    print("Conclusion: Par22 cannot interact with E3ub-insert105.")

    # 3. Calculate mass increase
    insertion_nucleotides = 105
    amino_acids = insertion_nucleotides / 3
    # Average MW of an amino acid is ~110 Daltons (0.11 kDa)
    mass_increase_kda = amino_acids * 0.11
    print("\n3. Insertion Mass Calculation:")
    print(f" - An insertion of {insertion_nucleotides} nucleotides in the coding sequence corresponds to {int(amino_acids)} amino acids ({insertion_nucleotides} / 3).")
    print(f" - The estimated mass increase is {int(amino_acids)} amino acids * ~110 Da/aa = {int(amino_acids*110)} Da, which is approximately {mass_increase_kda:.1f} kDa.")
    print("Conclusion: The insertion increases the mass of E3ub by approximately 4.0 kDa.")

    # 4. Calculate offspring resistance percentage
    self_pollination_rate = 0.05
    cross_pollination_rate = 0.95
    # For self-pollination (WI x WI), resistant offspring are WI (1/2) and II (1/4). Total = 3/4
    prob_resistant_self = 3/4
    # For cross-pollination (WI x WW), resistant offspring is WI (1/2). Total = 1/2
    prob_resistant_cross = 1/2
    
    total_prob_resistant = (prob_resistant_self * self_pollination_rate) + \
                           (prob_resistant_cross * cross_pollination_rate)
    percentage_resistant = total_prob_resistant * 100

    print("\n4. Theoretical Offspring Resistance Calculation:")
    print("The parent plant is heterozygous for the resistance allele (Genotype: WI). Resistance is conferred by the 'I' allele.")
    print(f"There is a {self_pollination_rate*100}% chance of self-pollination and a {cross_pollination_rate*100}% chance of cross-pollination with the wild-type population (Genotype: WW).")
    print("\nThe final equation for the percentage of resistant offspring is:")
    print(f"Percentage = [ (Prob. Resistant from Selfing × Selfing Rate) + (Prob. Resistant from Crossing × Crossing Rate) ] × 100")
    print(f"Percentage = [ ({prob_resistant_self} × {self_pollination_rate}) + ({prob_resistant_cross} × {cross_pollination_rate}) ] × 100")
    print(f"Percentage = [ {prob_resistant_self * self_pollination_rate} + {prob_resistant_cross * cross_pollination_rate} ] × 100")
    print(f"Percentage = [ {total_prob_resistant} ] × 100")
    print(f"Percentage = {percentage_resistant:.2f}%")
    
    print("\n---")
    print("Final Conclusion:")
    print("Based on the analysis, the correct statement must include:")
    print(f"- A theoretical resistance rate of {percentage_resistant:.2f}%.")
    print("- Only E3ub-wt is an active ubiquitin ligase.")
    print("- Par22 cannot interact with E3ub-insert105.")
    print("- A mass increase of approximately 4.0 kDa.")
    print("This matches option J.")

solve_biology_problem()
<<<J>>>