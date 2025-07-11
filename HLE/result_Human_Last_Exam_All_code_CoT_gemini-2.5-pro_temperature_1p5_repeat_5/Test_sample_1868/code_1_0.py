import sys

def solve_fair_score():
    """
    Calculates the theoretical maximum FAIR compliance score R based on a derived model.
    """
    # Given parameters
    c = 0.95  # Identifier resolution consistency
    b = 3     # Semantic versioning branching factor

    print("This script calculates the theoretical maximum FAIR score (R) for a federated knowledge graph system.")
    print("A model is constructed based on the FAIR principles and the given parameters.\n")
    print(f"Input Parameters:")
    print(f"  - Identifier Resolution Consistency (c): {c}")
    print(f"  - Semantic Versioning Branching Factor (b): {b}\n")
    print("--------------------------------------------------")
    print("Step 1: Modeling each FAIR principle (on a 0-1 scale)")
    print("--------------------------------------------------")

    # Model for Findability (f)
    # The maximum findability is limited by the probability of successfully resolving an identifier.
    f = c
    print(f"Findability (f) is limited by identifier consistency 'c'.")
    print(f"  f = c = {f}\n")

    # Model for Accessibility (a)
    # Accessibility also depends on successful identifier resolution.
    a = c
    print(f"Accessibility (a) is also limited by 'c', as access requires a valid identifier.")
    print(f"  a = c = {a}\n")

    # Model for Interoperability (i)
    # A higher branching factor 'b' leads to more semantic divergence, reducing interoperability.
    # An inverse relationship is a reasonable model.
    i = 1 / b
    print(f"Interoperability (i) is modeled as inversely proportional to the branching factor 'b'.")
    print(f"  i = 1 / b = 1 / {b} = {i:.4f}\n")

    # Model for Reusability (r)
    # Reusability requires findable, accessible, and interoperable data. It is limited by both 'c' and 'b'.
    r = c / b
    print(f"Reusability (r) depends on both identifier consistency and semantic stability.")
    print(f"  r = c / b = {c} / {b} = {r:.4f}\n")

    print("--------------------------------------------------")
    print("Step 2: Calculating the final FAIR Score (R)")
    print("--------------------------------------------------")

    # The overall score R is the scaled average of the four metrics.
    # The formula is: R = 10 * (f + a + i + r) / 4
    R = 10 * (f + a + i + r) / 4

    print("The final score R is the average of the four scores, scaled to 10.\n")
    print("Final Equation:")
    # Using sys.stdout.write to prevent an extra newline and format the output equation clearly
    sys.stdout.write(f"  R = 10 * (f + a + i + r) / 4\n")
    sys.stdout.write(f"  R = 10 * ({f} + {a} + {i:.4f} + {r:.4f}) / 4\n")
    sys.stdout.write(f"  R = 10 * ({f + a + i + r:.4f}) / 4\n")
    sys.stdout.write(f"  R = {10 * (f + a + i + r):.4f} / 4\n")
    print("--------------------------------------------------")
    print(f"Theoretical Maximum FAIR Score (R) = {R:.4f}")
    print("--------------------------------------------------")

solve_fair_score()