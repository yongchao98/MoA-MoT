import sys

# Define a function to solve the problem
def calculate_fair_score():
    """
    Calculates the theoretical maximum FAIR compliance score (R) for a 
    federated knowledge graph system based on given parameters.
    """
    # Given parameters
    # c: Consistency level of the decentralized identifier resolution mechanism
    c = 0.95
    # b: Branching factor of semantic version control
    b = 3

    # --- Model for FAIR Components (on a scale of 0 to 1) ---

    # Findability (f) is limited by the identifier consistency.
    f_max = c

    # Accessibility (a) is assumed to be perfect in a 'best implementation'.
    a_max = 1.0

    # Interoperability (i) is inversely affected by the semantic branching factor.
    i_max = 1 / b

    # Reusability (r) is dependent on both identifier consistency (for provenance)
    # and semantic stability (inversely related to branching).
    r_max = c * (1 / b)

    # The overall FAIR score (R) is the average of the four components.
    R_max = (f_max + a_max + i_max + r_max) / 4

    # --- Output Results ---
    
    print("This script calculates the theoretical maximum FAIR score (R) based on a defined model.")
    print(f"Given parameters: Identifier Consistency (c) = {c}, Semantic Branching (b) = {b}\n")
    
    print("The final equation is R = (f + a + i + r) / 4, where:")
    print(f"f (Findability) = c = {f_max:.4f}")
    print(f"a (Accessibility) = 1.0 (ideal) = {a_max:.4f}")
    print(f"i (Interoperability) = 1 / b = {i_max:.4f}")
    print(f"r (Reusability) = c * (1 / b) = {r_max:.4f}\n")

    # The final print statement showing the full calculation
    print("Plugging the values into the equation:")
    print(f"R = ({f_max:.2f} + {a_max:.1f} + {1/b:.2f} + {c * (1/b):.2f}) / 4")
    
    final_sum = f_max + a_max + i_max + r_max
    print(f"R = {final_sum:.2f} / 4")
    print(f"R = {R_max}")

    # The final answer in the required format for the platform.
    # We append '<<<...>>>' to the last line of the output stream.
    sys.stdout.write(f"\n<<<{R_max}>>>")

# Execute the function
calculate_fair_score()