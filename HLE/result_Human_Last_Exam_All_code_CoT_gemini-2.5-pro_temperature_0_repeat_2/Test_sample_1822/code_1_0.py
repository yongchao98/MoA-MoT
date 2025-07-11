import sys

def calculate_evidence_ratio():
    """
    Calculates the minimum ratio of evidence (E_r) required to confirm a true
    violation of fairness and not a statistical artifact.
    """
    # Given parameters from the problem description
    c = 0.85  # Confidence of the discriminatory pattern P in subgraph G
    s = 0.12  # Support of the pattern P in subgraph G
    c_prime = 0.78  # Confidence of the anti-correlated pattern in the 2-hop neighborhood
    beta = 0.23  # Bottleneck coefficient of the knowledge graph

    # The evidence for a violation is modeled as the product of its local confidence and support.
    # This captures both the strength and prevalence of the observed pattern.
    evidence_for_violation = c * s

    # The evidence for a statistical artifact is the counter-evidence from the neighborhood (c'),
    # modulated by the bottleneck coefficient (beta). A stronger bottleneck (lower beta)
    # reduces the relevance of the neighborhood evidence.
    evidence_for_artifact = c_prime * beta

    # The required evidence ratio (E_r) is the ratio of the evidence for the violation
    # to the evidence for it being a statistical artifact.
    if evidence_for_artifact == 0:
        # Avoid division by zero, though unlikely with the given inputs.
        # In this case, any positive evidence for violation would be infinitely stronger.
        if evidence_for_violation > 0:
            e_r = float('inf')
        else:
            e_r = 0.0 # Or undefined, but 0 is a safe numerical value.
    else:
        e_r = evidence_for_violation / evidence_for_artifact

    # Output the explanation and the calculation as requested.
    print("The minimum ratio of evidence (E_r) is calculated by comparing the evidence for the violation against the evidence for it being a statistical artifact.")
    print("\nFormula: E_r = (c * s) / (c' * β)")
    print("\nCalculation with the given values:")
    # The final equation with each number outputted
    print(f"E_r = ({c} * {s}) / ({c_prime} * {beta})")
    print(f"E_r = {evidence_for_violation:.4f} / {evidence_for_artifact:.4f}")
    print(f"E_r = {e_r:.4f}")

    # The final answer is returned to be captured by the wrapper
    return e_r

if __name__ == '__main__':
    final_answer = calculate_evidence_ratio()
    # The problem asks for the final answer in a specific format at the end of the response.
    # This is handled outside the script, but we can print it here for clarity.
    # For example: print(f"\n<<<{final_answer:.4f}>>>")