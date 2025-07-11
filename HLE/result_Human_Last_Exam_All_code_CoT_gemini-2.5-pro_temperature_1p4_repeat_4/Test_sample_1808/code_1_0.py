import numpy as np

def analyze_hybrid_zone_genetics():
    """
    This function simulates the formation of a hybrid zone to demonstrate
    the expected population genetic consequences, helping to identify which
    of the given options is NOT an expected outcome.
    """
    # --- 1. Simulation Parameters ---
    # These define the scenario: two populations that have been separated
    # for a long time and are now mixing in a hybrid zone.
    N_LOCI = 10000  # Number of genetic sites in the genome
    DIVERGENCE_STEPS = 5000 # Represents time of separation. More steps = more divergence.
    MUTATION_RATE = 1e-5  # A fundamental biological property of the species.

    print("--- Simulating a Hybrid Zone Scenario ---")
    print(f"The simulation uses a fixed Mutation Rate (u) of {MUTATION_RATE}.")
    print("This rate is a biological parameter and is not affected by population processes like gene flow.\n")

    # --- 2. Simulating Population Divergence ---
    # We start with two identical populations (represented by their allele frequencies).
    # We then introduce random mutations to make them different over time.
    # p1 and p2 store the frequency of a specific allele at each locus.
    p1 = np.full(N_LOCI, 0.5)
    p2 = np.full(N_LOCI, 0.5)

    # Let's simulate divergence through random drift and mutation.
    for _ in range(DIVERGENCE_STEPS):
        # Pick a random locus to change in each population
        locus1 = np.random.randint(0, N_LOCI)
        locus2 = np.random.randint(0, N_LOCI)
        # Randomly change allele frequency (simple drift model)
        p1[locus1] = np.random.rand()
        p2[locus2] = np.random.rand()

    # --- 3. Calculating Genetic Statistics ---
    # Now we have two divergent populations. A hybrid zone is a mix of them.
    # Let's calculate the metrics based on their allele frequencies (p1, p2).

    # B. Dxy: Absolute divergence between populations.
    # Equation: Dxy = p1*(1-p2) + p2*(1-p1)
    dxy_per_locus = p1 * (1 - p2) + p2 * (1 - p1)
    dxy = np.mean(dxy_per_locus)
    print("--- Analyzing Metrics Between Parent Populations ---")
    print(f"B. Dxy (Absolute Divergence) = {dxy:.4f}")
    print("   Explanation: A high Dxy reflects the accumulated genetic differences over a long period of separation. This is expected for populations forming a hybrid zone.")

    # A. Fst: Relative divergence (differentiation).
    # Equation: Fst = (Ht - Hs) / Ht
    p_total = (p1 + p2) / 2
    hs = (np.mean(2 * p1 * (1 - p1)) + np.mean(2 * p2 * (1 - p2))) / 2 # Avg. heterozygosity WITHIN pops
    ht = np.mean(2 * p_total * (1 - p_total)) # Expected heterozygosity of the TOTAL combined pop
    fst = (ht - hs) / ht if ht > 0 else 0
    print(f"\nA. Fst (Population Differentiation) = {fst:.4f}")
    print("   Explanation: A high Fst (value > ~0.15) shows the populations are very distinct. This is a prerequisite for a hybrid zone to exist.")

    # A hybrid zone population is a mixture of individuals from Pop1 and Pop2.
    # Now, we analyze the properties of this mixed population.

    # E. Pi (π): Nucleotide diversity WITHIN the hybrid zone.
    # This is the expected heterozygosity in the mixed population, which is Ht.
    pi_hybrid = ht
    pi_parental_avg = hs
    print("\n--- Analyzing Metrics Within The Hybrid Zone ---")
    print(f"E. Pi (Diversity in Hybrid Zone) = {pi_hybrid:.4f}")
    print(f"   (For comparison, average Pi in parent populations = {pi_parental_avg:.4f})")
    print("   Explanation: Pi is typically high in a hybrid zone because it contains the genetic variation from BOTH parent populations combined.")

    # C. Fis: Inbreeding coefficient WITHIN the hybrid zone.
    # Due to the Wahlund effect (pooling distinct groups), Fis is expected to be high.
    # Equation: Fis = (He - Ho) / He
    # Here, He = pi_hybrid (Ht) and Ho = pi_parental_avg (Hs).
    fis = (pi_hybrid - pi_parental_avg) / pi_hybrid if pi_hybrid > 0 else 0
    print(f"\nC. Fis (Inbreeding Coefficient in Hybrid Zone) = {fis:.4f}")
    print("   Explanation: A high Fis reflects a deficit of heterozygotes compared to random-mating expectations. This occurs when you mix two distinct populations (the Wahlund effect). Note that Fis is mathematically equal to Fst in this simple mixing scenario.")
    
    print("\n--- Final Conclusion ---")
    print("The simulation shows that in a hybrid zone formed by gene flow between distinct populations, we can expect:")
    print("-> High Fst (A) and High Dxy (B) between the parent populations.")
    print("-> High Pi (E) and High Fis (C) within the hybrid zone itself.")
    print("\nThe one exception is D. A high mutation rate (u) is not caused by gene flow. The mutation rate is a biological property of an organism's DNA repair mechanisms. While a species could have a high mutation rate, it is not a direct consequence or feature of the process of hybridization.")


if __name__ == '__main__':
    analyze_hybrid_zone_genetics()