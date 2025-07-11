def simulate_hybrid_zone_fis():
    """
    Calculates Fis in a simple, idealized hybrid zone to illustrate the effect of gene flow.

    Assumptions:
    - Parental Population 1 is fixed for allele 'A'.
    - Parental Population 2 is fixed for allele 'a'.
    - The hybrid zone population is formed by random mating between these two populations,
      resulting in F1 offspring.
    """
    print("Step 1: Define the parental populations and resulting hybrid offspring.")
    print("Parental Population 1 has genotype 'AA'.")
    print("Parental Population 2 has genotype 'aa'.")
    print("Offspring in the hybrid zone (F1) will all have the genotype 'Aa'.\n")

    # Step 2: Calculate allele frequencies in the hybrid zone
    # In a population of all 'Aa' individuals, the frequency of 'A' (p) is 0.5
    # and the frequency of 'a' (q) is 0.5.
    p = 0.5
    q = 0.5
    print(f"Step 2: Calculate allele frequencies in the hybrid population.")
    print(f"Frequency of allele 'A' (p) = {p}")
    print(f"Frequency of allele 'a' (q) = {q}\n")

    # Step 3: Calculate Expected Heterozygosity (He) under Hardy-Weinberg Equilibrium
    # He = 2 * p * q
    He = 2 * p * q
    print(f"Step 3: Calculate Expected Heterozygosity (He) using the formula He = 2 * p * q.")
    print(f"He = 2 * {p} * {q} = {He}\n")

    # Step 4: Determine Observed Heterozygosity (Ho)
    # Since all individuals are 'Aa', the frequency of heterozygotes is 1.0.
    Ho = 1.0
    print(f"Step 4: Determine Observed Heterozygosity (Ho).")
    print(f"In this idealized F1 hybrid population, all individuals are heterozygotes ('Aa').")
    print(f"Therefore, Ho = {Ho}\n")

    # Step 5: Calculate Fis
    # Fis = (He - Ho) / He
    Fis = (He - Ho) / He
    print(f"Step 5: Calculate Fis using the formula Fis = (He - Ho) / He.")
    print(f"The final equation is: Fis = ({He} - {Ho}) / {He}")
    print(f"Fis = {Fis}\n")

    print("Conclusion: Gene flow between distinct populations leads to an excess of heterozygotes,")
    print(f"resulting in a negative Fis of {Fis}. Therefore, a high positive Fis cannot occur.")

# Run the simulation
simulate_hybrid_zone_fis()