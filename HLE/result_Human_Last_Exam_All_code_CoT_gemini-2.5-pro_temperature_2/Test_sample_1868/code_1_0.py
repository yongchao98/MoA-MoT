import math

def calculate_fair_score():
    """
    Calculates the theoretical maximum FAIR compliance score (R) for a
    federated knowledge graph system based on given parameters.
    """
    # Given parameters
    c = 0.95  # Consistency level of decentralized identifiers
    b = 3     # Branching factor of semantic version control

    # --- Step 1: Model each FAIR component ---

    # Findability (f) is modeled as equal to the identifier consistency c.
    f = c

    # Accessibility (a) is assumed to be perfect (1.0) in a best-case scenario.
    a = 1.0

    # Interoperability (i) is modeled as inversely related to the complexity
    # introduced by the branching factor b.
    i = 1 / math.sqrt(b)

    # Reusability (r) is dependent on both findability and interoperability.
    # It is modeled as the product of f and i.
    r = f * i

    # --- Step 2: Calculate the overall FAIR score R ---

    # The overall score R is the arithmetic mean of the four components.
    R = (f + a + i + r) / 4

    # --- Step 3: Output the results ---

    print("--- FAIR Score Calculation ---")
    print(f"Given parameters: c = {c}, b = {b}\n")
    print("Modeled FAIR Components:")
    print(f"f (Findability)    = c                 = {f:.4f}")
    print(f"a (Accessibility)  = 1.0               = {a:.4f}")
    print(f"i (Interoperability) = 1 / sqrt({b})         = {i:.4f}")
    print(f"r (Reusability)    = f * i             = {r:.4f}\n")

    print("The final equation for the FAIR score R is:")
    print(f"R = (f + a + i + r) / 4")
    print(f"R = ({f:.4f} + {a:.4f} + {i:.4f} + {r:.4f}) / 4")

    print(f"\nThe theoretical maximum value of R is: {R:.4f}")

if __name__ == "__main__":
    calculate_fair_score()
    # Recalculating for the final answer tag with high precision
    # R = (0.95 + 1.0 + (1/math.sqrt(3)) + 0.95*(1/math.sqrt(3))) / 4
    # R = (1.95 + 1.95 * (1/math.sqrt(3))) / 4
    # R = 0.76895819...
    final_answer = 0.7690
    print(f"\n<<<{final_answer:.4f}>>>")
