def analyze_hybrid_zone_genetics():
    """
    Analyzes which population genetic metric is not an expected feature
    of a hybrid zone with ongoing gene flow.
    """

    # Step 1: Define Parent Populations and Prerequisites
    # A hybrid zone exists because two populations are genetically distinct.
    print("--- 1. Prerequisites for a Hybrid Zone ---")
    print("A hybrid zone forms between two populations that have diverged over time.")
    print("This implies:\n")

    # Analyze Fst
    print("A. High Fst (Fixation Index):")
    print("   - Fst measures genetic differentiation between populations.")
    print("   - High differentiation (High Fst) is the reason a hybrid zone can exist.")
    print("   - CONCLUSION: High Fst CAN occur.\n")

    # Analyze Dxy
    print("B. High Dxy (Absolute Divergence):")
    print("   - Dxy measures the average number of DNA sequence differences between populations, reflecting divergence time.")
    print("   - Like Fst, a high Dxy indicates a long period of separation and is a precondition for a distinct hybrid zone.")
    print("   - CONCLUSION: High Dxy CAN occur.\n")

    # Step 2: Analyze Consequences of Gene Flow Within the Hybrid Zone
    # We model the hybrid zone as a mix of individuals from two source populations.
    hybrid_zone_sample = ["Parent1_individual", "Parent2_individual", "F1_Hybrid_individual"]
    print("--- 2. Consequences within the Hybrid Zone ---")
    print(f"The zone itself contains a mix of individuals, such as: {hybrid_zone_sample}\n")

    # Analyze Pi
    print("E. High Pi (Nucleotide Diversity):")
    print("   - Pi measures the average genetic diversity *within* a population.")
    print("   - By mixing two divergent gene pools (from Parent 1 and Parent 2), the overall genetic diversity within the zone becomes high.")
    print("   - CONCLUSION: High Pi CAN occur.\n")

    # Analyze Fis
    print("C. High Fis (Inbreeding Coefficient):")
    print("   - Fis measures deviations from random mating within a population.")
    print("   - When distinct populations are pooled into one sample (the 'Wahlund Effect'), there is a statistical deficit of heterozygotes compared to what would be expected if it were a single, randomly-mating population.")
    print("   - A high positive Fis is a classic signature of a hybrid zone.")
    print("   - CONCLUSION: High Fis CAN occur.\n")

    # Step 3: Analyze the Independent Factor
    print("--- 3. Analyzing an Independent Evolutionary Force ---")

    # Analyze u
    print("D. High u (mutation rate):")
    print("   - 'u' is the fundamental rate at which new mutations arise in an organism's DNA due to, for example, replication errors.")
    print("   - Gene flow is the movement of *existing* genes (alleles) between populations.")
    print("   - There is no direct mechanism by which the process of inter-population mating alters the fundamental biochemical rate of mutation.")
    print("   - Gene flow and mutation are distinct evolutionary forces.")
    print("   - CONCLUSION: High 'u' CANNOT occur as a direct result of gene flow.\n")

    # Step 4: Final Summary "Equation"
    print("--- 4. Summary of Plausibility ---")
    print("Final Analysis Equation:")
    print("High Fst [A]: Possible (Prerequisite)")
    print("High Dxy [B]: Possible (Prerequisite)")
    print("High Fis [C]: Possible (Consequence)")
    print("High Pi  [E]: Possible (Consequence)")
    print("High u   [D]: Not a direct consequence of gene flow")


if __name__ == '__main__':
    analyze_hybrid_zone_genetics()
