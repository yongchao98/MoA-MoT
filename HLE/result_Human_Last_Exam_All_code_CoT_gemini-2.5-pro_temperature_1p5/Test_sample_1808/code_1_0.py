import numpy as np

def analyze_hybrid_zone():
    """
    Simulates a hybrid zone to demonstrate population genetic effects.
    This script shows that high Fst, Fis, and Pi can occur in a hybrid zone scenario,
    while the mutation rate (u) is not a direct consequence of gene flow.
    """

    # --- Step 1: Define Parental Populations ---
    # We model a single gene with two alleles, "A" and "a".
    # Population 1 is nearly fixed for allele "A".
    # Population 2 is nearly fixed for allele "a".
    # p is the frequency of allele "A", q is the frequency of allele "a".
    p1 = 0.99  # Freq(A) in Population 1
    q1 = 1 - p1 # Freq(a) in Population 1

    p2 = 0.01  # Freq(A) in Population 2
    q2 = 1 - p2 # Freq(a) in Population 2

    print("--- Analysis of a Hybrid Zone ---")
    print(f"Parental Population 1: Freq(A)={p1}, Freq(a)={q1:.2f}")
    print(f"Parental Population 2: Freq(A)={p2}, Freq(a)={q2:.2f}")
    print("-" * 40)

    # --- Step 2: Analyze Fst (Option A) ---
    # Fst measures differentiation. High Fst is expected between the source populations.
    # Fst = (Ht - Hs) / Ht
    # Ht: Expected heterozygosity in the total (combined) population
    # Hs: Average expected heterozygosity across the two source populations
    p_total = (p1 + p2) / 2
    q_total = (q1 + q2) / 2
    
    # Expected heterozygosities (2pq)
    Ht = 2 * p_total * q_total
    Hs1 = 2 * p1 * q1
    Hs2 = 2 * p2 * q2
    Hs = (Hs1 + Hs2) / 2
    
    # Avoid division by zero if Ht is 0
    Fst = (Ht - Hs) / Ht if Ht > 0 else 0

    print("A. Can High Fst occur? YES.")
    print("Fst measures differentiation between the source populations.")
    print(f"Fst = ({Ht:.4f} - {Hs:.4f}) / {Ht:.4f} = {Fst:.4f}")
    print("Result: Fst is very high, as expected for distinct populations forming a hybrid zone.")
    print("-" * 40)

    # --- Step 3: Analyze Fis (Option C) ---
    # Fis measures heterozygote deficit within the hybrid zone due to mixing (Wahlund effect).
    # When populations first mix, observed heterozygosity is just the average of the parental populations (Hs).
    # The expected heterozygosity for the combined gene pool is Ht.
    # Fis = (H_expected - H_observed) / H_expected
    H_observed_in_hybrid_zone = Hs
    H_expected_in_hybrid_zone = Ht
    Fis = (H_expected_in_hybrid_zone - H_observed_in_hybrid_zone) / H_expected_in_hybrid_zone if H_expected_in_hybrid_zone > 0 else 0
    
    print("C. Can High Fis occur? YES.")
    print("Fis measures heterozygote deficit within the new hybrid population.")
    print(f"Fis = ({H_expected_in_hybrid_zone:.4f} - {H_observed_in_hybrid_zone:.4f}) / {H_expected_in_hybrid_zone:.4f} = {Fis:.4f}")
    print("Result: Fis is very high, indicating a strong heterozygote deficit (Wahlund effect).")
    print("-" * 40)

    # --- Step 4: Analyze Pi (Option E) ---
    # Pi (nucleotide diversity) is equivalent to expected heterozygosity.
    # We compare Pi in the hybrid zone to the parental populations.
    pi_hybrid_zone = Ht # Diversity in the combined gene pool
    pi_parent1 = Hs1
    pi_parent2 = Hs2
    
    print("E. Can High Pi occur? YES.")
    print("Pi (genetic diversity) can be higher in the hybrid zone than in either parent.")
    print(f"Pi in Parent 1: {pi_parent1:.4f}")
    print(f"Pi in Parent 2: {pi_parent2:.4f}")
    print(f"Pi in Hybrid Zone: {pi_hybrid_zone:.4f}")
    print("Result: Pi is significantly higher in the hybrid zone due to the mixing of alleles.")
    print("-" * 40)

    # --- Step 5: Final Conclusion about all options ---
    print("Summary:")
    print("A. High Fst: Occurs between source populations.")
    print("B. High Dxy: Occurs, reflects historical divergence between source populations.")
    print("C. High Fis: Occurs within the hybrid zone due to the Wahlund effect.")
    print("E. High Pi: Occurs within the hybrid zone due to the mixing of diverse alleles.")
    print("\nD. High u (mutation rate): CANNOT occur as a direct result of gene flow.")
    print("Mutation rate (u) is a biochemical property, not a population-level consequence of allele mixing.")

if __name__ == '__main__':
    analyze_hybrid_zone()