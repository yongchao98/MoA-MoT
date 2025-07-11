import sys

def solve_fair_score():
    """
    Calculates the theoretical maximum FAIR compliance score (R) based on
    identifier consistency (c) and semantic branching factor (b).
    """
    # Given parameters
    c = 0.95  # Identifier resolution consistency level
    b = 3     # Semantic version control branching factor

    # --- Step 1: Explain the model ---
    print("This script calculates the theoretical maximum FAIR score (R) for a federated knowledge graph.")
    print("The model assumes R is the average of its four components: Findability (f), Accessibility (a), Interoperability (i), and Reusability (r).")
    print("\nInitial Values:")
    print(f"Identifier Consistency (c) = {c}")
    print(f"Branching Factor (b) = {b}")

    # --- Step 2: Model the impact of constraints ---
    # Findability and Accessibility are limited by identifier consistency.
    f_max = c
    a_max = c

    # Interoperability and Reusability are limited by branching complexity.
    i_max = 1 / b
    r_max = 1 / b
    
    print("\nModeling the Maximum Score for Each FAIR Component:")
    print(f"Max Findability (f_max) = c = {f_max:.4f}")
    print(f"Max Accessibility (a_max) = c = {a_max:.4f}")
    print(f"Max Interoperability (i_max) = 1/b = {i_max:.4f}")
    print(f"Max Reusability (r_max) = 1/b = {r_max:.4f}")

    # --- Step 3: Calculate the final theoretical maximum R ---
    R_max = (f_max + a_max + i_max + r_max) / 4.0

    print("\nCalculating the Theoretical Maximum FAIR Score (R_max):")
    # As requested, printing the equation with the final numbers
    print(f"R_max = (f_max + a_max + i_max + r_max) / 4")
    print(f"R_max = ({f_max:.4f} + {a_max:.4f} + {i_max:.4f} + {r_max:.4f}) / 4")
    
    final_result_str = f"{R_max:.4f}"
    print(f"\nThe final theoretical maximum value of R is: {final_result_str}")
    
    # Appending the final answer in the specified format to the output stream.
    # We use sys.stdout.write to avoid adding an extra newline.
    sys.stdout.write(f"\n<<<{final_result_str}>>>")

# Execute the solution
solve_fair_score()