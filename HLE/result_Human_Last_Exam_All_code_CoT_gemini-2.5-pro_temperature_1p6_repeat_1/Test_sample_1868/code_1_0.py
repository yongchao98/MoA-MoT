import sys

def calculate_fair_score():
    """
    Calculates the theoretical maximum FAIR compliance score (R) for a federated knowledge graph.

    This model assumes:
    1. R is the arithmetic mean of its four components (f, a, i, r).
    2. Findability (f) and Accessibility (a) are limited by the identifier resolution consistency (c).
    3. Interoperability (i) and Reusability (r) are limited by both consistency (c) and penalized
       by the semantic versioning branching factor (b), modeled as a 1/b degradation.
    """
    # Given parameters
    c = 0.95  # Identifier resolution consistency level
    b = 3     # Semantic version control branching factor

    # --- Step 1: Model the theoretical maximums for each FAIR component ---
    # Findability is limited by identifier consistency
    f = c
    # Accessibility is also limited by identifier consistency
    a = c
    # Interoperability is limited by consistency and semantic branching
    i = c / b
    # Reusability is limited by consistency and semantic branching
    r = c / b

    # --- Step 2: Calculate the final FAIR score R ---
    R = (f + a + i + r) / 4

    # --- Step 3: Print the breakdown and the final result ---
    print("This script calculates the theoretical maximum FAIR score (R) for a federated knowledge graph system.")
    print("\nGiven Parameters:")
    print(f"Identifier Consistency (c): {c}")
    print(f"Branching Factor (b): {b}")
    
    print("\nFAIR Component Score Calculation:")
    print(f"Findability (f) = c = {f:.4f}")
    print(f"Accessibility (a) = c = {a:.4f}")
    print(f"Interoperability (i) = c / b = {c:.2f} / {b} = {i:.4f}")
    print(f"Reusability (r) = c / b = {c:.2f} / {b} = {r:.4f}")

    print("\nFinal FAIR Score (R) Calculation:")
    # Per instructions, outputting each number in the final equation
    print(f"R = (f + a + i + r) / 4")
    print(f"R = ({f:.4f} + {a:.4f} + {i:.4f} + {r:.4f}) / 4")
    
    # Python 3.8+ f-string can have self-documenting expressions
    if sys.version_info >= (3, 8):
        print(f"R = {(f + a + i + r):.4f} / 4 = {R:.4f}")
    else:
        sum_components = f + a + i + r
        print(f"R = {sum_components:.4f} / 4 = {R:.4f}")


if __name__ == "__main__":
    calculate_fair_score()
    # The final answer as requested
    c = 0.95
    b = 3
    R = (c + c + c/b + c/b) / 4
    # To conform to the specified output format, print the value inside <<<>>>
    print(f"\n<<<{R:.4f}>>>")