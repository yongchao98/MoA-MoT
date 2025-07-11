import math

def calculate_fair_score():
    """
    Calculates the theoretical maximum FAIR compliance score (R) for a
    federated knowledge graph system based on given parameters.
    """
    # Given parameters
    c = 0.95  # Consistency level of decentralized identifier resolution
    b = 3     # Branching factor of semantic version control

    # Step 1: Model and calculate each FAIR component score based on the plan.
    # Findability (f) is limited by the consistency level.
    f = c

    # Accessibility (a) is assumed to be perfect (1.0) in the best-case scenario.
    a = 1.0

    # Interoperability (i) is influenced by consistency (c) and penalized by branching (b).
    i = c / b

    # Reusability (r) is also influenced by consistency (c) and penalized by branching (b).
    r = c / b

    # Step 2: Calculate the final FAIR score R as the average of the components.
    R = (f + a + i + r) / 4

    # Step 3: Print the breakdown of the calculation and the final result.
    print("Calculating the theoretical maximum FAIR score (R):")
    print("-" * 50)
    print(f"Given parameters:")
    print(f"  Consistency Level (c): {c}")
    print(f"  Branching Factor (b): {b}\n")
    print("Component Scores:")
    print(f"  f (Findability)   = c                  = {f}")
    print(f"  a (Accessibility) = 1.0 (ideal)        = {a}")
    print(f"  i (Interoperability)= c / b = {c} / {b} = {i}")
    print(f"  r (Reusability)     = c / b = {c} / {b} = {r}\n")

    print("Final Calculation using R = (f + a + i + r) / 4:")
    # The final equation with each number explicitly shown
    print(f"R = ({f} + {a} + {i} + {r}) / 4")
    print(f"\nThe theoretical maximum value of R is: {R}")

calculate_fair_score()
print("\n<<<0.6458333333333333>>>")