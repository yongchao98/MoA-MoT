import sys

# Redirect print to a string to control output format.
# This part is for internal logic and will not be in the final code block.
# The final code block will use standard print.

# This function helps in creating a clear, final output.
def evaluate_experiments():
    """
    Analyzes the provided honeybee experiment data to determine the correct conclusion.
    """
    # --- Data from the problem description ---
    baseline_mortality = 10  # Mortality rate of non-infected bees in %

    # General mortality for Fungus A (most pollens, excluding sunflower's effect)
    fungus_a_general_mortality = 35

    # General mortality for Fungus B (consistent across all pollens)
    fungus_b_general_mortality = 20

    # General mortality for Fungus C (consistent across all pollens)
    fungus_c_general_mortality = 10

    # --- Step 1: Analyze Pathogenicity of Fungus A and B ---
    print("--- Analysis of Fungus A and B ---")
    print(f"The baseline mortality for non-infected bees is {baseline_mortality}%.")

    # Check Fungus A
    print(f"The general mortality rate for bees infected with Fungus A is {fungus_a_general_mortality}%.")
    if fungus_a_general_mortality > baseline_mortality:
        print(f"Since {fungus_a_general_mortality} > {baseline_mortality}, Fungus A is a pathogen.")
    else:
        print(f"Since {fungus_a_general_mortality} <= {baseline_mortality}, Fungus A is not a pathogen.")

    # Check Fungus B
    print(f"The mortality rate for bees infected with Fungus B is {fungus_b_general_mortality}%.")
    if fungus_b_general_mortality > baseline_mortality:
        print(f"Since {fungus_b_general_mortality} > {baseline_mortality}, Fungus B is a pathogen.")
    else:
        print(f"Since {fungus_b_general_mortality} <= {baseline_mortality}, Fungus B is not a pathogen.")
    print("\nConclusion 1: Both Fungus A and Fungus B are pathogens as they increase mortality.\n")


    # --- Step 2: Analyze Fungus C ---
    print("--- Analysis of Fungus C ---")
    print(f"The mortality rate for bees infected with Fungus C is {fungus_c_general_mortality}%.")
    if fungus_c_general_mortality > baseline_mortality:
        print(f"Since {fungus_c_general_mortality} > {baseline_mortality}, Fungus C is a pathogen.")
    else:
        print(f"Since the mortality rate of {fungus_c_general_mortality}% is not higher than the baseline of {baseline_mortality}%, Fungus C is not a pathogen.")
        print("Furthermore, Experiment 5 shows that Fungus C infection often increases egg production, suggesting a non-harmful relationship.")
        print("Therefore, Fungus C can be classified as a commensal (or even symbiotic) organism.")

    print("\nConclusion 2: Fungus C is a commensal.\n")

    # --- Step 3: Final Conclusion based on Analysis ---
    print("--- Final Evaluation ---")
    print("Based on the analysis:")
    print("1. Fungus A and B are pathogens.")
    print("2. Fungus C is a commensal.")
    print("This directly corresponds to answer choice I.")


# Run the analysis
evaluate_experiments()