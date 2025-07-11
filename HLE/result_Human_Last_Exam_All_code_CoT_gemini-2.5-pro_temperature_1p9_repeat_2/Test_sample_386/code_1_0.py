import sys

def analyze_culture_conditions(fbs_percent, antibiotic_percent):
    """
    Analyzes the suitability of medium components for culturing corneal fibroblasts.

    Args:
        fbs_percent (int): The percentage of Fetal Bovine Serum (FBS).
        antibiotic_percent (int): The percentage of antibiotic solution.
    """
    # Step 1: Explain the biological context of the chosen option.
    print("The task is to identify the correct method for establishing a fibroblast cell line from the cornea.")
    print("Option C describes a biologically sound and standard procedure.\n")
    print("This involves:")
    print("- Isolating stromal cells after removing the epithelium and endothelium.")
    print("- Culturing these cells in a medium that promotes growth and adhesion.")
    print("- Observing their differentiation into myofibroblasts, a common in vitro event.\n")

    # Step 2: Validate and explain the specific numerical parameters.
    # The prompt asks to output each number in the final equation.
    print(f"To achieve this, the culture medium composition, or 'equation', is crucial:")
    print(f"Final Medium = Basal Medium + {fbs_percent}% FBS + {antibiotic_percent}% Antibiotic\n")
    
    # Step 3: Justify why these percentages are correct.
    is_fbs_correct = (fbs_percent == 10)
    is_antibiotic_correct = (antibiotic_percent == 1)

    if is_fbs_correct:
        print(f"Validation of Serum: {fbs_percent}% FBS is a standard concentration.")
        print(" -> It provides the necessary growth factors for fibroblast proliferation and their differentiation into myofibroblasts.\n")
    else:
        print(f"Validation of Serum: {fbs_percent}% FBS is not the typical concentration cited.", file=sys.stderr)

    if is_antibiotic_correct:
        print(f"Validation of Antibiotic: {antibiotic_percent}% is a standard concentration.")
        print(" -> It effectively prevents bacterial contamination without being toxic to the cells.\n")
    else:
        print(f"Validation of Antibiotic: {antibiotic_percent}% is not a typical concentration.", file=sys.stderr)

    if is_fbs_correct and is_antibiotic_correct:
        print("Conclusion: The parameters in Option C are correct for the described procedure.")


# Parameters from the correct option (C)
fbs_percentage_c = 10
antibiotic_percentage_c = 1

# Run the analysis
analyze_culture_conditions(fbs_percentage_c, antibiotic_percentage_c)
