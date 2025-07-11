# peptide_synthesis_advisor.py

def suggest_synthesis_technique():
    """
    Analyzes the synthesis challenge for a long peptide with an unnatural amino acid
    and recommends the most suitable technique.
    """

    # --- Problem Definition ---
    peptide_length = 100
    unnatural_amino_acid = "Azido Phenylalanine (pAzF)"
    sequence_snippet = "M...[KAVCLXVIGATR[...]A"
    
    # --- Introduction ---
    print("--- Peptide Synthesis Strategy Advisor ---")
    print(f"\nTask: Recommend a synthesis method for a {peptide_length}aa peptide.")
    print(f"Special requirement: Site-specific incorporation of an unnatural amino acid, '{unnatural_amino_acid}'.")

    # --- Analysis of Techniques ---
    print("\n--- Analysis of Potential Techniques ---")

    # 1. Solid-Phase Peptide Synthesis (SPPS)
    print("\n[1] Solid-Phase Peptide Synthesis (SPPS)")
    print("    - Description: Step-by-step chemical addition of amino acids on a solid support.")
    print("    - Suitability: While it can incorporate unnatural amino acids, it is not ideal for long peptides.")
    print(f"    - Challenge: For a {peptide_length}aa peptide, the cumulative yield becomes extremely low.")
    print("      (e.g., 99% efficiency per step results in 0.99^99 ≈ 37% theoretical max yield, with practical yields being much lower).")
    print("    - Verdict: Impractical due to low yield and significant purification challenges.")

    # 2. Native Chemical Ligation (NCL)
    print("\n[2] Native Chemical Ligation (NCL)")
    print("    - Description: A chemical reaction that joins two smaller peptide fragments.")
    print("    - Suitability: Excellent for building long peptides. The 100aa peptide could be made by joining two ~50aa fragments.")
    print("    - Strategy: One fragment with the UAA can be made via SPPS, and the other can be made recombinantly.")
    print(f"    - Note: This is a very strong option, especially since the sequence snippet '{sequence_snippet}' contains a Cysteine ('C'), which is required for NCL.")
    print("    - Verdict: A highly viable and flexible strategy.")

    # 3. Recombinant Expression using an Expanded Genetic Code
    print("\n[3] Recombinant Expression using an Expanded Genetic Code")
    print("    - Description: Uses an engineered organism (like E. coli) to produce the entire peptide.")
    print("    - How it works:")
    print("      a) An orthogonal tRNA/synthetase pair is introduced into the host cell.")
    print(f"      b) This pair specifically recognizes the UAA ({unnatural_amino_acid}) and a unique codon (e.g., the UAG amber stop codon).")
    print("      c) The UAA is supplied in the growth media and incorporated by the ribosome during translation.")
    print("    - Suitability: This technique is designed precisely for producing large proteins with site-specifically incorporated UAAs.")
    print("    - Verdict: Often the most direct and high-yielding method for this specific problem.")

    # --- Final Recommendation ---
    print("\n--- Recommendation ---")
    print("While NCL is a powerful strategy, the single most helpful and efficient technique for producing a full-length 100aa peptide with a UAA like Azido Phenylalanine is:")
    print("\n********************************************************")
    print("*  Recombinant expression using an expanded genetic code *")
    print("********************************************************")


# Execute the advisor function
if __name__ == '__main__':
    suggest_synthesis_technique()