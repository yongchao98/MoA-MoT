def analyze_fmr1_gene():
    """
    Analyzes and explains the categories of FMR1 gene CGG repeats
    to determine the likely molecular abnormality in the patient.
    """

    # Define the CGG repeat count thresholds for the FMR1 gene
    NORMAL_MAX = 44
    INTERMEDIATE_MAX = 54
    PREMUTATION_MAX = 200

    # The minimums are derived from the maximum of the previous category
    INTERMEDIATE_MIN = NORMAL_MAX + 1
    PREMUTATION_MIN = INTERMEDIATE_MAX + 1
    FULL_MUTATION_MIN = PREMUTATION_MAX + 1

    print("Analysis of FMR1 Gene Trinucleotide Repeats")
    print("---------------------------------------------")
    print(f"The FMR1 gene's function is related to the number of its CGG repeats.")
    print("\nThere are four main categories:")

    # Print each category with its numerical range
    print(f"1. Normal Range: 5-{NORMAL_MAX} repeats. No associated health problems.")
    print(f"2. Intermediate Range: {INTERMEDIATE_MIN}-{INTERMEDIATE_MAX} repeats. Generally no symptoms, but may have a slightly increased risk for some issues.")
    print(f"3. Premutation Range: {PREMUTATION_MIN}-{PREMUTATION_MAX} repeats. This is the range relevant to the patient.")
    print(f"4. Full Mutation Range: >{PREMUTATION_MAX} repeats. Causes Fragile X Syndrome.")

    print("\nPatient Correlation:")
    print("The patient presents with ovarian dysgenesis, amenorrhea, and infertility with a normal 46,XX karyotype.")
    print("This clinical picture is characteristic of Fragile X-associated Primary Ovarian Insufficiency (FXPOI).")
    print(f"\nConclusion:")
    print(f"FXPOI is caused by a premutation in the FMR1 gene. Therefore, the likely molecular abnormality is a CGG repeat expansion into the premutation range.")
    print(f"The patient likely has between {PREMUTATION_MIN} and {PREMUTATION_MAX} CGG repeats in one of her FMR1 genes.")

# Execute the analysis
analyze_fmr1_gene()
