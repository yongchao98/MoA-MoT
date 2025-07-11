import sys

def calculate_fair_score():
    """
    Calculates the theoretical maximum FAIR compliance score (R) for a
    federated knowledge graph system based on given parameters.
    """
    # Given parameters for the federated knowledge graph system F
    c = 0.95  # Consistency level of the decentralized identifier resolution
    b = 3     # Branching factor of semantic version control

    # --- Step 1: Model the FAIR metrics based on the parameters ---
    # Findability (f) is limited by the consistency of identifier resolution.
    f = c
    # Accessibility (a) is also limited by the ability to resolve identifiers.
    a = c
    # Interoperability (i) is inversely affected by semantic branching.
    # More branches lead to fragmentation, reducing interoperability.
    i = 1 / b
    # Reusability (r) is affected by both identifier consistency (for provenance)
    # and semantic stability.
    r = c / b

    # --- Step 2: Calculate the overall FAIR score R as the mean of the components ---
    R = (f + a + i + r) / 4

    # --- Step 3: Print the detailed calculation process ---
    print("Calculating the theoretical maximum FAIR score (R)")
    print("-" * 50)
    print(f"Input Parameters:")
    print(f"  - Identifier Resolution Consistency (c): {c}")
    print(f"  - Semantic Versioning Branching Factor (b): {b}\n")

    print("Modeled FAIR Component Scores:")
    print(f"  - Findability (f) = c = {f:.4f}")
    print(f"  - Accessibility (a) = c = {a:.4f}")
    print(f"  - Interoperability (i) = 1/b = {i:.4f}")
    print(f"  - Reusability (r) = c/b = {r:.4f}\n")

    print("Final Score Calculation (R):")
    # As requested, the final equation with each number is printed below
    print(f"  R = (f + a + i + r) / 4")
    print(f"  R = ({f:.4f} + {a:.4f} + {i:.4f} + {r:.4f}) / 4")
    print(f"  R = {f + a + i + r:.4f} / 4")
    print(f"  R = {R:.4f}")
    
    # Use stdout for the final answer to be captured
    sys.stdout.write(f"\n<<<{R:.4f}>>>\n")

if __name__ == '__main__':
    calculate_fair_score()