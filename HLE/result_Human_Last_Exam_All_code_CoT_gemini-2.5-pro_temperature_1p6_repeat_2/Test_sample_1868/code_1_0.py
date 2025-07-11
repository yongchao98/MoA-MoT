def calculate_fair_score():
    """
    Calculates the theoretical maximum FAIR compliance score (R) for a federated knowledge graph.

    This calculation is based on the following model:
    - R = (f + a + i + r) / 4
    - Findability (f) is determined by identifier consistency (c).
    - Accessibility (a) is assumed to be perfect (1.0) under best practices.
    - Interoperability (i) and Reusability (r) are impacted by both identifier
      consistency (c) and semantic ambiguity from the versioning branching factor (b).
    """
    # Given parameters
    c = 0.95  # Consistency level of identifier resolution
    b = 3     # Branching factor of semantic version control

    # Modeling each FAIR component based on the plan
    f = c
    a = 1.0   # Assumed perfect accessibility in the best-case scenario
    i = c / b
    r = c / b

    # Calculate the final FAIR compliance score R
    R = (f + a + i + r) / 4

    # --- Output ---
    print("Step 1: Define the model for each FAIR component.")
    print(f"Findability (f) = c = {f}")
    print(f"Accessibility (a) = {a} (assumed best practice)")
    print(f"Interoperability (i) = c / b = {c} / {b} = {i:.4f}")
    print(f"Reusability (r) = c / b = {c} / {b} = {r:.4f}")
    print("\nStep 2: Calculate the final FAIR score R using R = (f + a + i + r) / 4.")
    # The final equation with all numbers shown explicitly as requested
    print("\nFinal Equation:")
    print(f"R = ({f} + {a} + {i:.4f} + {r:.4f}) / 4")
    
    # Calculate and print the final result
    print(f"R = {R:.4f}")


if __name__ == "__main__":
    calculate_fair_score()
