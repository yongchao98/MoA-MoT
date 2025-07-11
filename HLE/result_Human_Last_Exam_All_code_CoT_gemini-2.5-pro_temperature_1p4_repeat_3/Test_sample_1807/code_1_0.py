import sys

def analyze_genetic_differentiation():
    """
    Analyzes and explains genetic differentiation between sexes
    by simulating an autosomal and a sex-linked marker.
    """

    print("Analyzing potential causes for genetic differentiation between males and females...")
    print("-" * 70)

    # --- Case 1: Autosomal Marker ---
    print("Case 1: An Autosomal Marker (e.g., on Chromosome 1)")
    print("In a randomly mating population, autosomal markers are inherited independently of sex.")
    print("Therefore, allele frequencies are expected to be equal in males and females.")

    # Example frequencies for an allele 'A'
    autosomal_freq_male = 0.6
    autosomal_freq_female = 0.6

    print(f"   - Frequency of allele 'A' in males:   {autosomal_freq_male}")
    print(f"   - Frequency of allele 'A' in females: {autosomal_freq_female}")
    # Fst is calculated as (Ht - Hs) / Ht. If allele frequencies are identical, Fst = 0.
    print("   - Resulting Fst: 0 (No differentiation)")
    print("-" * 70)

    # --- Case 2: Sex-Linked Marker (Y-Chromosome) ---
    print("Case 2: A Sex-Linked Marker (on the Y chromosome in an XY system)")
    print("The Y chromosome is only present in males. Therefore, any marker unique to it will only be found in males.")

    # Allele 'P' (for present) on the Y chromosome.
    # It exists in 100% of males and 0% of females.
    y_linked_freq_male = 1.0
    y_linked_freq_female = 0.0

    print(f"   - Frequency of Y-linked marker in males:   {y_linked_freq_male}")
    print(f"   - Frequency of Y-linked marker in females: {y_linked_freq_female}")
    # With frequencies of 1 and 0, the populations are completely distinct at this locus.
    print("   - Resulting Fst: 1 (Maximum possible differentiation)")
    print("-" * 70)

    # --- Conclusion ---
    print("\nConclusion:")
    print("The simulation shows that while autosomal markers exhibit no differentiation between sexes,")
    print("markers on sex chromosomes (like the Y chromosome) show extreme differentiation.")
    print("This is a direct consequence of the genetic makeup of sex determination systems (XY or ZW).")
    print("Therefore, finding that some markers exhibit pronounced differentiation is a classic signature")
    print("of markers located on sex chromosomes.")

# It seems that the script execution is not correctly redirecting the stdout.
# A temporary solution is to write the output to a file and then print the file content.
original_stdout = sys.stdout
with open("output.txt", "w") as f:
    sys.stdout = f
    analyze_genetic_differentiation()
sys.stdout = original_stdout

with open("output.txt", "r") as f:
    print(f.read())
