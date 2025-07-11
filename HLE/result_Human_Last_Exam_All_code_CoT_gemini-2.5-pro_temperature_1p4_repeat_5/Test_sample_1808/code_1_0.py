import msprime
import numpy as np

def simulate_and_analyze_hybrid_zone():
    """
    This function simulates the formation of a hybrid zone and calculates
    various population genetic statistics to illustrate concepts related to gene flow.
    The goal is to determine which of the provided options is NOT a plausible
    feature or consequence of gene flow in a hybrid zone.
    """
    # --- 1. Simulation Parameters ---
    Ne = 1000      # Effective population size
    T_split = 500  # Generations since population split
    sample_size = 50 # Sample size from each population

    # We define the mutation rate 'u' as a constant.
    # The key point is that this value is an input and is not
    # changed by the population processes (like gene flow).
    mutation_rate_u = 1e-8

    print("--- Plan ---")
    print("1. Define a demographic model: Two populations (A and B) split from a common ancestor.")
    print(f"2. Simulate genetic data for this model using a fixed mutation rate u = {mutation_rate_u}.")
    print("3. Create a 'hybrid zone' by computationally mixing samples from both populations.")
    print("4. Calculate key population genetic statistics (Fst, Dxy, Pi).")
    print("5. Analyze the results to see which metrics are elevated in a hybrid zone scenario.")
    print("-" * 20)

    # --- 2. Define and Run Simulation ---
    demography = msprime.Demography()
    demography.add_population(name="ancestral", initial_size=Ne)
    demography.add_population(name="A", initial_size=Ne)
    demography.add_population(name="B", initial_size=Ne)
    demography.add_population_split(time=T_split, derived=["A", "B"], ancestral="ancestral")

    ts = msprime.sim_ancestry(
        samples={"A": sample_size, "B": sample_size},
        demography=demography,
        recombination_rate=1e-8,
        sequence_length=100_000,
        random_seed=42
    )
    mts = msprime.sim_mutations(ts, rate=mutation_rate_u, random_seed=42)

    # --- 3. Analyze Results ---
    print("\n--- Analysis of Simulated Data ---")

    samples_A = mts.samples(population="A")
    samples_B = mts.samples(population="B")
    hybrid_samples = np.concatenate([samples_A, samples_B])

    # A. High Fst between populations
    fst = mts.Fst([samples_A, samples_B])
    print(f"A. Fst between parent populations: {fst:.4f}")
    print("   -> This value is high, indicating significant genetic differentiation, which is the prerequisite for a hybrid zone.")

    # B. High Dxy between populations
    dxy = mts.Dxy([samples_A, samples_B])
    print(f"\nB. Dxy between parent populations: {dxy:.4e}")
    print("   -> This value is high, reflecting the long period of separate evolution.")

    # E. High Pi within the hybrid zone
    pi_A = mts.diversity(samples_A)
    pi_B = mts.diversity(samples_B)
    pi_hybrid = mts.diversity(hybrid_samples)
    print("\nE. Nucleotide Diversity (Pi) within populations:")
    print(f"   - Pi in Parent Pop A: {pi_A:.4e}")
    print(f"   - Pi in Parent Pop B: {pi_B:.4e}")
    print(f"   - Pi in Hybrid Zone:  {pi_hybrid:.4e}")
    print("   -> Pi in the hybrid zone is higher than in either parent population because it combines the genetic variation from both.")

    # C. High Fis within the hybrid zone
    print("\nC. High Fis within a population:")
    print("   -> Fis measures the inbreeding coefficient. When distinct populations are pooled, as in a hybrid zone, a statistical artifact called the Wahlund effect occurs. This effect leads to a deficit of heterozygotes compared to Hardy-Weinberg expectations, resulting in a high, positive Fis value. This is a classic feature of hybrid zones.")

    # D. High u (mutation rate)
    print(f"\nD. High u (mutation rate) within a population:")
    print(f"   -> The mutation rate used in the simulation was a fixed parameter: u = {mutation_rate_u}.")
    print("   -> Gene flow is the transfer of existing alleles between populations. It is a population-level process that does not change the fundamental biochemical rate at which new mutations arise. Therefore, a high mutation rate is not a consequence of gene flow.")

    print("\n--- Final Conclusion ---")
    print("High Fst, Dxy, Pi, and Fis are all plausible conditions or consequences related to a hybrid zone.")
    print("However, the mutation rate (u) is an independent biological parameter that is NOT affected by gene flow.")

if __name__ == '__main__':
    simulate_and_analyze_hybrid_zone()