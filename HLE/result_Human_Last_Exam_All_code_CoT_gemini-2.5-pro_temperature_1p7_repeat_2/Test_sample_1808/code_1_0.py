import sys
import io

# Helper to redirect print to a string for analysis if needed,
# but for this case, direct printing is fine.

def analyze_hybrid_zone_genetics():
    """
    Simulates the effect of gene flow on the Fis statistic in a hybrid zone
    to determine which of the multiple-choice options cannot occur.
    """
    print("This analysis demonstrates the effect of gene flow on population genetic statistics.")
    print("=" * 70)

    # --- Initial State: A Hybrid Zone Forms (Wahlund Effect) ---
    # Two populations, P1 (all AA) and P2 (all aa), meet.
    # We sample the zone as a single population made of 50% P1 and 50% P2 individuals
    # before they have had a chance to interbreed.

    print("Step 1: Analyzing the hybrid zone *before* widespread gene flow occurs.")
    print("The population is a mix of two distinct groups, causing a Wahlund Effect.\n")

    # Overall allele frequencies (p and q) in the mixed population
    p_A = 0.5  # Frequency of allele 'A'
    q_a = 0.5  # Frequency of allele 'a'

    # Expected Heterozygosity (He) based on Hardy-Weinberg Equilibrium
    # He = 2 * p * q
    He = 2 * p_A * q_a

    # Observed Heterozygosity (Ho) in the non-interbred mix is 0,
    # because there are only AA and aa individuals, and no Aa individuals.
    Ho_initial = 0.0

    # Calculate Fis = (He - Ho) / He
    # A positive Fis indicates a deficit of heterozygotes.
    if He == 0:
        Fis_initial = 0
    else:
        Fis_initial = (He - Ho_initial) / He

    print("Initial State Calculation:")
    print(f"  - Allele 'A' frequency (p): {p_A}")
    print(f"  - Allele 'a' frequency (q): {q_a}")
    print(f"  - Expected Heterozygosity (He = 2*p*q): 2 * {p_A} * {q_a} = {He}")
    print(f"  - Observed Heterozygosity (Ho): {Ho_initial} (no interbreeding yet)")
    print(f"  - Initial Fis = (He - Ho) / He: ({He} - {Ho_initial}) / {He} = {Fis_initial:.2f}")
    print("\nResult: A high Fis of 1.00 is observed, indicating a complete lack of heterozygotes.")
    print("-" * 70)


    # --- After Gene Flow: Random Mating Occurs ---
    # Now, we simulate one generation of gene flow (i.e., random mating) in the zone.
    # The new generation's genotypes will follow HWE proportions.

    print("Step 2: Simulating one generation of gene flow (random mating).")
    print("This mating process creates hybrid individuals (heterozygotes).\n")

    # The new Observed Heterozygosity (Ho) will be equal to the
    # expected frequency of heterozygotes (2pq).
    Ho_after_mating = 2 * p_A * q_a

    # Recalculate Fis with the new Ho. He remains the same as allele frequencies are unchanged.
    if He == 0:
        Fis_new = 0
    else:
        Fis_new = (He - Ho_after_mating) / He

    print("Post-Gene-Flow Calculation:")
    print(f"  - Expected Heterozygosity (He) remains: {He}")
    print(f"  - New Observed Heterozygosity (Ho) is now equal to He: {Ho_after_mating}")
    print(f"  - New Fis = (He - Ho) / He: ({He} - {Ho_after_mating}) / {He} = {Fis_new:.2f}")
    print("\nResult: Fis drops to 0.00, indicating heterozygosity matches random mating expectations.")
    print("=" * 70)

    # --- Conclusion ---
    print("\nConclusion:")
    print("The process of gene flow involves mating between groups, which *creates* heterozygotes.")
    print("This simulation shows that gene flow actively reduces a high Fis towards zero.")
    print("Therefore, a sustained 'High Fis' (heterozygote deficit) is contradictory to the process of active gene flow.")
    print("\nReview of other options:")
    print(" A. High Fst: Plausible. The source populations remain distinct.")
    print(" B. High Dxy: Plausible. Reflects historical divergence between populations.")
    print(" D. High u (mutation): Plausible. Mutation rate is an independent process.")
    print(" E. High Pi (diversity): Plausible. Gene flow increases genetic diversity (Pi) within the zone.")

analyze_hybrid_zone_genetics()