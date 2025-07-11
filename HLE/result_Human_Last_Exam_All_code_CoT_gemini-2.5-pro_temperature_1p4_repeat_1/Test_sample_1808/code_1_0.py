import sys

def calculate_wahlund_effect(p1, p2):
    """
    Calculates and explains the Wahlund effect for a hybrid zone
    formed by two populations with different allele frequencies.

    Args:
      p1 (float): Frequency of allele 'A' in population 1.
      p2 (float): Frequency of allele 'A' in population 2.
    """
    print("This script demonstrates the Wahlund effect, which causes high Fis in a hybrid zone.\n")

    # Allele frequencies for q
    q1 = 1 - p1
    q2 = 1 - p2

    # In a 50/50 mix, the overall allele frequencies in the hybrid zone are the average
    p_total = (p1 + p2) / 2
    q_total = (q1 + q2) / 2

    # 1. Calculate Expected Heterozygosity (He) in the hybrid zone
    # This is what we'd expect if the zone was one large, randomly mating population.
    he_total = 2 * p_total * q_total

    # 2. Calculate Observed Heterozygosity (Ho) in the hybrid zone
    # This is the average of the heterozygosity from the two source populations.
    he1 = 2 * p1 * q1
    he2 = 2 * p2 * q2
    ho_total = (he1 + he2) / 2 # Assuming 50/50 mix

    # 3. Calculate Fis
    # Fis = (Expected_He - Observed_Ho) / Expected_He
    # A positive Fis indicates a deficit of heterozygotes.
    if he_total == 0:
        fis = 0
    else:
        fis = (he_total - ho_total) / he_total

    print("--- Setup ---")
    print(f"Population 1 allele 'A' frequency (p1): {p1}")
    print(f"Population 2 allele 'A' frequency (p2): {p2}")

    print("\n--- Calculation ---")
    print(f"Overall allele 'A' frequency in hybrid zone (p_total): {p_total:.2f}")
    print(f"Expected Heterozygosity (He) = 2 * p_total * q_total = 2 * {p_total:.2f} * {q_total:.2f} = {he_total:.2f}")
    print(f"Observed Heterozygosity (Ho) = (He_pop1 + He_pop2) / 2 = ({he1:.2f} + {he2:.2f}) / 2 = {ho_total:.2f}")

    print("\n--- Result ---")
    print("Fis is the inbreeding coefficient. A high positive value indicates a heterozygote deficit.")
    sys.stdout.write("Fis = (He - Ho) / He = ")
    # The prompt requests printing the equation with numbers
    print(f"({he_total:.2f} - {ho_total:.2f}) / {he_total:.2f} = {fis:.4f}")

    print("\nBecause Fis is high, this demonstrates that high Fis (Option C) CAN occur.")
    print("The only option that is not a direct consequence of gene flow in a hybrid zone is a change in the mutation rate (u).")


# Let's model two divergent populations
# Pop 1 is nearly fixed for allele 'A'
p1_freq = 0.9
# Pop 2 is nearly fixed for allele 'a'
p2_freq = 0.1

calculate_wahlund_effect(p1_freq, p2_freq)
<<<D>>>