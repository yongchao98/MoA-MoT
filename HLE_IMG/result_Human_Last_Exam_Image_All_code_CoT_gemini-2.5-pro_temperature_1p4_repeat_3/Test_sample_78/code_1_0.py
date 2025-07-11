def analyze_zebrafish_genotypes():
    """
    Analyzes a gel image from an RFLP experiment to determine zebrafish genotypes.

    The function is based on the interpretation of the provided gel image where
    a heterozygous cross results in three distinct banding patterns corresponding
    to the three possible genotypes.

    Genotype Interpretation based on the gel:
    - Homozygous Wild-Type (+/+): Has only the uncut PCR product (one high band).
    - Homozygous Mutant (m/m): Has only the cut PCR product (two low bands).
    - Heterozygote (+/m): Has both uncut and cut products (three bands).
    """

    # Based on visual inspection of the 18 sample lanes in the provided gel image:
    
    # Count of lanes with only the single, high-molecular-weight (uncut) band.
    # These are lanes 2 and 18.
    num_homozygous_wild_type = 2

    # Count of lanes with all three bands (one high, two low).
    # These are lanes 1, 6, 7, 8, 9, 10, 11, 12, 13, 14, 17.
    num_heterozygotes = 11

    # Count of lanes with only the two, lower-molecular-weight (cut) bands.
    # These are lanes 3, 4, 5, 15, 16.
    num_homozygous_mutants = 5
    
    # Print the result in the format A/B/C, where:
    # A = # of homozygous wild-type larvae
    # B = # of heterozygotes
    # C = # of homozygous mutants
    
    print("The counts for each genotype are:")
    print(f"Homozygous Wild-Type (A): {num_homozygous_wild_type}")
    print(f"Heterozygotes (B): {num_heterozygotes}")
    print(f"Homozygous Mutants (C): {num_homozygous_mutants}")
    print("\nFinal Answer in A/B/C format:")
    print(f"{num_homozygous_wild_type}/{num_heterozygotes}/{num_homozygous_mutants}")

analyze_zebrafish_genotypes()