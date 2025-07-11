import sys

def calculate_risk_change():
    """
    This script models the change in disease risk based on an HLA variant.
    """
    # The problem states the fold increase in antigen presentation probability.
    fold_increase = 1000

    # Let's assume a hypothetical, arbitrary base probability of presentation
    # for an individual without the variant. We'll use a small number.
    base_probability = 0.00005

    # Calculate the new, increased probability for an individual with the variant.
    new_probability = base_probability * fold_increase

    # The biological conclusion is that increased presentation of a disease-causing
    # self-antigen leads to a higher risk of autoimmune disease.
    print("An HLA variant that increases the presentation of a disease-causing self-antigen would likely increase a person's risk of developing the disease.")
    print("\n--- Calculation Model ---")
    print("This happens because the immune system is more likely to see and react to the self-antigen.")
    
    # As requested, we will print the final equation showing each number.
    print("\nFinal Equation:")
    # The f-string formats the floating-point numbers for clarity.
    print(f"New Probability = Base Probability * Fold Increase")
    print(f"{new_probability:.4f} = {base_probability} * {fold_increase}")

if __name__ == '__main__':
    calculate_risk_change()
