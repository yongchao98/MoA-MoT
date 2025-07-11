def analyze_zebrafish_genotypes():
    """
    Analyzes a virtual gel image to determine zebrafish genotypes
    based on a PCR-RFLP experiment.
    """
    print("Step 1: Understanding the RFLP genotyping logic.")
    print("The experiment uses PCR followed by a restriction digest with SfaNI.")
    print("A C-to-A mutation is present in the mutant allele.")
    print("We assume this mutation destroys a natural SfaNI site found in the wild-type allele.\n")

    print("Step 2: Predicting the banding patterns for each genotype.")
    print(" - Homozygous Wild-Type (+/+): The SfaNI site is present on both alleles. The PCR product is cut, resulting in two smaller bands.")
    print(" - Homozygous Mutant (-/-): The SfaNI site is absent from both alleles. The PCR product is not cut, resulting in one large band.")
    print(" - Heterozygote (+/-): One allele is cut and one is not. This results in three bands (the two small cut fragments and the one large uncut fragment).\n")

    print("Step 3: Counting the patterns on the provided gel image.")
    # These counts are based on careful visual inspection of the 15 sample lanes in the gel image.
    # Pattern A (two bands, representing +/+)
    # Pattern B (three bands, representing +/-)
    # Pattern C (one band, representing -/-)
    num_two_bands = 6  # Lanes 2, 3, 4, 5, 13, 14
    num_three_bands = 7 # Lanes 1, 7, 8, 9, 10, 11, 12
    num_one_band = 2   # Lanes 6, 15
    
    total_larvae = num_one_band + num_two_bands + num_three_bands
    print(f"Based on the gel image, I counted {total_larvae} individual larvae with the following patterns:")
    print(f" - Number of lanes with two bands (Homozygous Wild-Type): {num_two_bands}")
    print(f" - Number of lanes with three bands (Heterozygote): {num_three_bands}")
    print(f" - Number of lanes with one band (Homozygous Mutant): {num_one_band}\n")

    # Assigning counts based on the logic above
    A = num_two_bands
    B = num_three_bands
    C = num_one_band

    print("Step 4: Final genotype counts.")
    print(f"The number of homozygous wild-type larvae (A) is {A}.")
    print(f"The number of heterozygous larvae (B) is {B}.")
    print(f"The number of homozygous mutant larvae (C) is {C}.\n")

    print("Final Answer in 'A/B/C' format:")
    print(f"{A}/{B}/{C}")


analyze_zebrafish_genotypes()
<<<6/7/2>>>