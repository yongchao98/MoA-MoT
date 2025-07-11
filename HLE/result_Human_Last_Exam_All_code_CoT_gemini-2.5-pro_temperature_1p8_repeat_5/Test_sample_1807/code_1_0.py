def calculate_male_female_fst():
    """
    Calculates and explains the Fst between males and females for an X-linked gene.

    This demonstrates how sex chromosome systems can lead to high differentiation
    between the sexes at specific markers.
    """
    # Assume a gene on the X chromosome with two alleles, 'A' and 'a'.
    # Let's say due to selection or drift, the frequency of allele 'A' has
    # diverged between males and females.
    # p_females is the frequency of allele 'A' in females (who are XX)
    # p_males is the frequency of allele 'A' in males (who are XY)
    p_females = 0.8
    p_males = 0.3

    print(f"Scenario: An X-linked marker where the frequency of an allele is {p_females} in females and {p_males} in males.\n")

    # Fst can be calculated as Var(p) / (p_bar * (1 - p_bar))
    # Where p_bar is the mean allele frequency and Var(p) is the variance
    # of allele frequencies between the subpopulations (males and females).

    # For X-linked genes, females (XX) contribute 2/3 of the X chromosomes to the
    # gene pool and males (XY) contribute 1/3. We use these weights.
    weight_f = 2/3
    weight_m = 1/3

    # 1. Calculate the weighted average allele frequency (p_bar)
    p_bar = (p_females * weight_f) + (p_males * weight_m)
    print("Step 1: Calculate the weighted average allele frequency (p_bar)")
    print(f"p_bar = (p_females * 2/3) + (p_males * 1/3)")
    print(f"p_bar = ({p_females} * {weight_f:.3f}) + ({p_males} * {weight_m:.3f}) = {p_bar:.4f}\n")

    # 2. Calculate the weighted variance of allele frequencies (var_p)
    var_p = weight_f * (p_females - p_bar)**2 + weight_m * (p_males - p_bar)**2
    print("Step 2: Calculate the weighted variance of allele frequencies (var_p)")
    print("var_p = (2/3)*(p_females - p_bar)^2 + (1/3)*(p_males - p_bar)^2")
    print(f"var_p = {weight_f:.3f}*({p_females} - {p_bar:.4f})^2 + {weight_m:.3f}*({p_males} - {p_bar:.4f})^2 = {var_p:.4f}\n")

    # 3. Calculate the total expected heterozygosity (Ht)
    Ht = p_bar * (1 - p_bar)
    print("Step 3: Calculate the expected heterozygosity in the total population (Ht)")
    print("Ht = p_bar * (1 - p_bar)")
    print(f"Ht = {p_bar:.4f} * (1 - {p_bar:.4f}) = {Ht:.4f}\n")

    # 4. Calculate Fst
    # Fst is undefined if Ht is 0 (i.e., if the allele is fixed)
    if Ht == 0:
        fst = float('nan')
        print("Total heterozygosity is 0, Fst is undefined.")
    else:
        fst = var_p / Ht
        print("Step 4: Calculate Fst")
        print("Fst = var_p / Ht")
        print(f"Fst = {var_p:.4f} / {Ht:.4f} = {fst:.4f}\n")

    print(f"The resulting Fst value of {fst:.4f} indicates pronounced genetic differentiation between males and females for this marker.")

calculate_male_female_fst()