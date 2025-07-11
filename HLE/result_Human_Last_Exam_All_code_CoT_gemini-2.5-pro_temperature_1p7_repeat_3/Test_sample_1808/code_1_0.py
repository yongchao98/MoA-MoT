def calculate_fst(p1, p2):
    """Calculates Fst for two populations with given allele frequencies."""
    # He or Hs is the average expected heterozygosity within subpopulations
    # He = 2 * p * q
    hs1 = 2 * p1 * (1 - p1)
    hs2 = 2 * p2 * (1 - p2)
    hs_avg = (hs1 + hs2) / 2.0

    # Ht is the expected heterozygosity in the total, pooled population
    p_total = (p1 + p2) / 2.0
    ht = 2 * p_total * (1 - p_total)

    # Fst formula
    if ht == 0:
        return 0 # No variation in the total population
    fst = (ht - hs_avg) / ht
    return fst

def main():
    """
    Simulates gene flow between two populations to show its effect on Fst and Pi.
    """
    # --- Step 1: Initial State (No Gene Flow) ---
    # Two completely diverged populations at a single locus
    # Population 1 has only allele 'A'
    p1_initial = 1.0
    # Population 2 has only allele 'a'
    p2_initial = 0.0

    # Pi (related to Heterozygosity, 2pq) within each pop is 0, as they are fixed.
    pi1_initial = 2 * p1_initial * (1 - p1_initial)
    pi2_initial = 2 * p2_initial * (1 - p2_initial)
    fst_initial = calculate_fst(p1_initial, p2_initial)

    print("--- Initial State (Before Gene Flow) ---")
    print(f"Population 1 Allele 'A' Frequency: {p1_initial}")
    print(f"Population 2 Allele 'A' Frequency: {p2_initial}")
    print(f"Initial Diversity (Pi proxy) in Pop 1: {pi1_initial}")
    print(f"Initial Diversity (Pi proxy) in Pop 2: {pi2_initial}")
    print(f"Initial Fst between populations: {fst_initial:.2f}")
    print("Interpretation: Fst = 1.0 means complete differentiation with zero gene flow.")
    print("\n" + "="*50 + "\n")


    # --- Step 2: Gene Flow Event ---
    # Let's model a migration rate (m) of 20% between the populations
    migration_rate = 0.2
    p1_final = p1_initial * (1 - migration_rate) + p2_initial * migration_rate
    p2_final = p2_initial * (1 - migration_rate) + p1_initial * migration_rate


    # --- Step 3: Final State (After Gene Flow) ---
    # Pi (diversity) within each population is now non-zero.
    pi1_final = 2 * p1_final * (1 - p1_final)
    pi2_final = 2 * p2_final * (1 - p2_final)
    fst_final = calculate_fst(p1_final, p2_final)

    print("--- Final State (After Gene Flow) ---")
    print(f"Population 1 Allele 'A' Frequency: {p1_final:.2f}")
    print(f"Population 2 Allele 'A' Frequency: {p2_final:.2f}")
    print(f"Final Diversity (Pi proxy) in Pop 1: {pi1_final:.2f}")
    print(f"Final Diversity (Pi proxy) in Pop 2: {pi2_final:.2f}")
    print(f"Final Fst between populations: {fst_final:.2f}")
    print("Interpretation: Gene flow has introduced new alleles, increasing diversity (Pi)")
    print("within each population and has reduced differentiation, lowering Fst.")


if __name__ == "__main__":
    main()