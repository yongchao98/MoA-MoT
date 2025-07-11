import sys

def calculate_max_fair_score():
    """
    Calculates the theoretical maximum FAIR compliance score (R) for a federated
    knowledge graph system based on given system parameters.
    """

    # --- Problem Parameters ---
    # c: Consistency level of the decentralized identifier resolution mechanism.
    c = 0.95
    # b: Branching factor of semantic version control for each knowledge graph.
    b = 3

    # --- Plan Explanation ---
    print("This script calculates the theoretical maximum FAIR score (R) based on a quantitative model.")
    print("The model assumes the following relationships:")
    print("1. Findability (f) and Accessibility (a) are capped by identifier consistency 'c'.")
    print("2. Interoperability (i) and Reusability (r) are capped by 'c' and penalized by the semantic branching factor 'b' (as c/b).")
    print("3. The final score R is the average of the four components, scaled to 10.")
    print("-" * 60)

    # --- Step 1: Model the maximum achievable scores for each FAIR principle (0-1 scale) ---
    f_max = c
    a_max = c
    i_max = c / b
    r_max = c / b
    
    print("Given Parameters:")
    print(f"  Identifier Consistency (c) = {c}")
    print(f"  Semantic Branching Factor (b) = {b}")
    print("\nMaximum Component Scores (0-1 Scale):")
    print(f"  Max Findability (f)     = c      = {f_max:.4f}")
    print(f"  Max Accessibility (a)   = c      = {a_max:.4f}")
    print(f"  Max Interoperability (i)  = c / b  = {i_max:.4f}")
    print(f"  Max Reusability (r)     = c / b  = {r_max:.4f}")
    print("-" * 60)
    
    # --- Step 2: Model and Calculate the final score R (0-10 scale) ---
    # The final score R is the scaled average of the four component scores.
    # R = 10 * (f + a + i + r) / 4
    R_max = 10 * (f_max + a_max + i_max + r_max) / 4
    
    # --- Step 3: Output the final equation and result ---
    print("The final score (R_max) is calculated using the formula:")
    print("R_max = 10 * (f_max + a_max + i_max + r_max) / 4\n")
    
    # Per the instructions, we output each number in the final equation.
    print("Substituting the component scores into the equation:")
    # We use round() for a cleaner display in the equation line, but use the full precision value for the final result.
    print(f"R_max = 10 * ({round(f_max, 4)} + {round(a_max, 4)} + {round(i_max, 4)} + {round(r_max, 4)}) / 4")

    # The exact value is 19/3
    print(f"\nTheoretical Maximum FAIR Score (R_max) = {R_max:.4f}")

# Execute the function
calculate_max_fair_score()
