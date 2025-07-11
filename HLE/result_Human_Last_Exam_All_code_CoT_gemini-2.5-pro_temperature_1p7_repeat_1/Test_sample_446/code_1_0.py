import sys

def analyze_integrin_binding_peptides():
    """
    Analyzes a list of RGD-containing peptides to determine which is most likely
    to bind an integrin receptor with high affinity based on known biological data.
    """

    # Dictionary containing the peptides and an explanation of their binding likelihood.
    # The explanation is based on established biochemical knowledge.
    peptide_info = {
        'A. RGDMAA': "Contains the core RGD motif, but the flanking 'MAA' is not a known high-affinity enhancer.",
        'B. RGDSPSS': "Strong candidate. This sequence is a close analog of the binding site in fibronectin (a major integrin ligand). The 'SP' (Serine-Proline) sequence following RGD is known to be critical for inducing a conformation that binds with high affinity to integrins like α5β1.",
        'C. RGDLTTP': "Contains the core RGD motif, but the flanking 'LTTP' is not a well-characterized high-affinity sequence.",
        'D. RGDQVSK': "Contains the core RGD motif, but the flanking 'QVSK' does not correspond to a known high-affinity binding enhancer.",
        'E. RGDARGG': "Contains the core RGD motif, but the flanking 'ARGG' is not characteristic of high-affinity integrin binding sites."
    }

    best_candidate_key = None
    # We identify the best candidate based on keywords indicating a known, potent biological motif.
    for key, description in peptide_info.items():
        if "fibronectin" in description or "critical" in description:
            best_candidate_key = key
            break

    print("Analysis of Peptides for Integrin Binding:")
    print("-" * 45)
    for key, description in peptide_info.items():
        print(f"Peptide {key}: {description}")
    print("-" * 45)
    
    if best_candidate_key:
        print("\nConclusion:")
        print(f"The peptide most likely to show strong binding in an in vitro assay is from choice {best_candidate_key.split('.')[0]}.")
        print("This is because its sequence closely mimics a well-studied, high-affinity binding site found in the natural integrin ligand, fibronectin.")
    else:
        print("\nCould not determine the best candidate based on the provided information.")


# Execute the analysis
analyze_integrin_binding_peptides()