import math

def calculate_evidence_ratio():
    """
    Calculates the minimum ratio of evidence (E_r) required to confirm a fairness violation.

    The model synthesizes local evidence (confidence 'c' and support 's'),
    contradictory neighborhood evidence (confidence 'c_prime'), and the graph's
    structural bottleneck coefficient ('beta').

    E_r represents the factor by which the current evidence for a violation
    must be strengthened to overcome the counter-evidence.
    """
    # Step 1: Define the given variables
    c = 0.85  # Confidence of discriminatory pattern in subgraph G
    s = 0.12  # Support of the pattern in G
    c_prime = 0.78  # Confidence of anti-correlated pattern in the 2-hop neighborhood
    beta = 0.23  # Bottleneck coefficient of the knowledge graph

    # Step 2: Calculate the components of the formula
    # The local evidence strength is discounted by the bottleneck coefficient.
    local_evidence_strength = c * s * beta

    # The required ratio E_r is the counter-evidence divided by the discounted local evidence.
    # E_r = c_prime / (c * s * beta)
    e_r = c_prime / local_evidence_strength

    # Step 3: Print the step-by-step calculation
    print("The formula for the minimum evidence ratio (E_r) is:")
    print("E_r = c' / (c * s * β)\n")

    print("Substituting the given values:")
    # The f-string formats the numbers into the equation string.
    # Using 'g' for general format to avoid unnecessary trailing zeros.
    print(f"E_r = {c_prime} / ({c} * {s} * {beta})")

    print("\nCalculating the denominator (the discounted local evidence strength):")
    print(f"E_r = {c_prime} / {c * s * beta}")

    print("\nFinal calculation:")
    print(f"E_r = {e_r}")

    # Return the final value for the 'answer' block
    return e_r

# Execute the function and capture the final answer
final_answer = calculate_evidence_ratio()

# Final answer block as requested by the user format
print(f"\n<<<{final_answer}>>>")