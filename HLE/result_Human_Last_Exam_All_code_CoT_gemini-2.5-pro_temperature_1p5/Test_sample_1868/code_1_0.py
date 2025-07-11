import sys

def calculate_fair_score():
    """
    Calculates the theoretical maximum FAIR compliance score (R) for a
    federated knowledge graph system based on given constraints.
    """
    # Given parameters from the problem description
    c = 0.95  # Consistency level of identifier resolution
    b = 3     # Branching factor of semantic version control

    # --- Step 1: Model the impact of constraints on FAIR components ---
    # Findability (f) and Accessibility (a) are limited by the consistency 'c'.
    # In the best case, their scores cannot exceed the reliability of the
    # underlying identifier system.
    f_max = c
    a_max = c

    # Interoperability (i) and Reusability (r) are limited by the ambiguity
    # introduced by the branching factor 'b'. The score is modeled as being
    # inversely proportional to the number of branches.
    i_max = 1 / b
    r_max = 1 / b

    # --- Step 2: Calculate the overall FAIR score R as the average ---
    # The total score R is the average of the four component scores.
    R_max = (f_max + a_max + i_max + r_max) / 4

    # --- Step 3: Print the detailed explanation and results ---
    print("Calculating the theoretical maximum FAIR compliance score (R)...")
    print("-" * 60)
    print(f"Given System Parameters:")
    print(f"  - Identifier Resolution Consistency (c): {c}")
    print(f"  - Semantic Versioning Branching Factor (b): {b}\n")

    print("Modeling the Maximum Scores for each FAIR Component:")
    print("  - Findability and Accessibility are capped by identifier consistency 'c'.")
    print("  - Interoperability and Reusability are penalized by branching factor 'b'.\n")

    print("Component Score Calculation:")
    print(f"  - Maximum Findability (f_max)    = c       = {f_max:.4f}")
    print(f"  - Maximum Accessibility (a_max)   = c       = {a_max:.4f}")
    print(f"  - Maximum Interoperability (i_max) = 1/b = 1/{b} = {i_max:.4f}")
    print(f"  - Maximum Reusability (r_max)      = 1/b = 1/{b} = {r_max:.4f}\n")

    print("Final FAIR Score (R) Calculation:")
    print("  R_max = (f_max + a_max + i_max + r_max) / 4")
    print(f"  R_max = ({f_max:.4f} + {a_max:.4f} + {i_max:.4f} + {r_max:.4f}) / 4")
    print(f"  R_max = {f_max + a_max + i_max + r_max:.4f} / 4")
    print(f"  R_max = {R_max:.4f}\n")

    # Final answer rounded to three decimal places
    final_answer = round(R_max, 3)
    
    # Writing the final answer to stdout in the required format
    # Redirecting original stdout to stderr to hide this print from normal output
    original_stdout = sys.stdout
    sys.stdout = sys.stderr
    print(f"<<<{final_answer}>>>")
    sys.stdout = original_stdout


if __name__ == "__main__":
    calculate_fair_score()