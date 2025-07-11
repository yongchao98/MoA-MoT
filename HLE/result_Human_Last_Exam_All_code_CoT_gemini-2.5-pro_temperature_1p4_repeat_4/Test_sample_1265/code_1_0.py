import collections

def analyze_allele_frequency(depth, pileup_string):
    """
    Analyzes a pileup string to calculate allele frequencies.

    Args:
        depth (int): The sequencing depth.
        pileup_string (str): A string representing the bases from each read at a specific position.
    """
    print(f"--- Analyzing at {depth}X Depth ---")
    print(f"Pileup data ({len(pileup_string)} reads): {pileup_string}")
    
    # Count the occurrences of each base (A, C, G, T)
    counts = collections.Counter(pileup_string)
    
    # Sort the counts for consistent output
    sorted_counts = sorted(counts.items())
    
    print("Allele Counts:")
    # The "equation" is the final count of each base
    equation_parts = []
    for base, count in sorted_counts:
        frequency = (count / depth) * 100
        print(f"  Base {base}: Count = {count}, Frequency = {frequency:.2f}%")
        equation_parts.append(f"{base}:{count}")

    print("\nFinal Count Equation:")
    print(" + ".join(equation_parts) + f" = {depth} total reads")
    print("-" * 25 + "\n")


def main():
    """
    Main function to run the demonstration.
    """
    # At 40X depth, we expect our 5% resistance allele ('T') to appear 2 times (0.05 * 40).
    # We'll add one random sequencing error ('G'). The rest is the wild-type ('A').
    pileup_40x = "A" * 37 + "T" * 2 + "G" * 1
    
    # At 80X depth, we expect our 5% resistance allele ('T') to appear 4 times (0.05 * 80).
    # We'll add one random sequencing error ('G'). The signal of the true allele is now stronger.
    pileup_80x = "A" * 75 + "T" * 4 + "G" * 1

    print("Demonstrating the effect of increased sequencing depth on identifying low-frequency alleles.")
    print("Wild-type allele: 'A'")
    print("True low-frequency resistance allele: 'T' (at 5%)")
    print("Random sequencing error: 'G'\n")

    # Analyze the low-depth scenario
    analyze_allele_frequency(40, pileup_40x)

    # Analyze the high-depth scenario
    analyze_allele_frequency(80, pileup_80x)

    print("Conclusion:")
    print("At 40X depth, the resistance allele 'T' (5.00%) and the error 'G' (2.50%) have similar low frequencies, making them hard to distinguish.")
    print("At 80X depth, the 'T' signal (5.00%) is much stronger relative to the 'G' error (1.25%), providing higher confidence that 'T' is a true allele.")

if __name__ == "__main__":
    main()
