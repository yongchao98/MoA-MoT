def calculate_fair_score():
    """
    Calculates the theoretical maximum FAIR compliance score (R) for a
    federated knowledge graph system based on given parameters.
    """
    # Given parameters
    c = 0.95  # Identifier resolution consistency level
    b = 3     # Semantic version control branching factor

    # --- Modeling the FAIR components ---

    # Maximum Findability (f) is limited by the identifier resolution consistency.
    max_f = c

    # Maximum Accessibility (a) is assumed to be perfect (1.0) in the best case.
    max_a = 1.0

    # Maximum Interoperability (i) is inversely affected by the branching factor,
    # as more branches create ambiguity.
    max_i = 1 / b

    # Maximum Reusability (r) is also inversely affected by the branching factor,
    # as it complicates provenance tracking.
    max_r = 1 / b

    # --- Calculating the theoretical maximum R score ---

    # The overall FAIR score R is the average of the four components.
    max_R = (max_f + max_a + max_i + max_r) / 4

    # --- Output the results ---
    
    # Print the equation with the values plugged in
    print("The calculation for the theoretical maximum FAIR score (R) is:")
    print(f"R = (f + a + i + r) / 4")
    print(f"R = ({max_f} + {max_a} + {max_i:.4f} + {max_r:.4f}) / 4")
    
    # Print the final result
    print(f"\nTheoretical Maximum R = {max_R:.4f}")

if __name__ == "__main__":
    calculate_fair_score()