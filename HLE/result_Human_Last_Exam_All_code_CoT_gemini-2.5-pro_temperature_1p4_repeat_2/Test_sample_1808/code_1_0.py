import math

def demonstrate_fst_vs_gene_flow():
    """
    This function demonstrates the inverse relationship between gene flow (m)
    and population differentiation (Fst).
    """
    print("This script will show why high Fst cannot occur with high gene flow.")
    print("We use the island model formula: Fst = 1 / (1 + 4 * N * m)")
    print("where N is effective population size and m is the migration rate (gene flow).\n")

    # Assume a constant effective population size
    N = 1000
    print(f"Let's assume an effective population size (N) of {N}.\n")

    # Test different rates of gene flow
    gene_flow_rates = [0, 0.0001, 0.001, 0.01, 0.1]

    for m in gene_flow_rates:
        # Calculate the denominator for clarity
        denominator = 1 + (4 * N * m)
        # Calculate Fst
        fst = 1 / denominator

        print(f"--- Calculating for a gene flow rate (m) of: {m} ---")
        print(f"Equation: Fst = 1 / (1 + 4 * N * m)")
        print(f"Calculation: Fst = 1 / (1 + 4 * {N} * {m})")
        print(f"Calculation: Fst = 1 / ({denominator})")
        print(f"Result: Fst = {fst:.4f}")

        if m == 0:
            print("Interpretation: With ZERO gene flow, Fst is at its maximum (1.0). Populations are completely differentiated.")
        elif fst > 0.25:
            print("Interpretation: This is a high Fst, indicating very limited gene flow.")
        else:
            print("Interpretation: This is a low Fst, resulting from significant gene flow.")
        print("")

    print("="*60)
    print("Conclusion: As the rate of gene flow (m) increases, the value of Fst decreases.")
    print("Therefore, high gene flow is incompatible with high Fst. When gene flow occurs, it actively reduces Fst.")
    print("This means high Fst between populations CANNOT occur when there is significant gene flow.")

if __name__ == '__main__':
    demonstrate_fst_vs_gene_flow()