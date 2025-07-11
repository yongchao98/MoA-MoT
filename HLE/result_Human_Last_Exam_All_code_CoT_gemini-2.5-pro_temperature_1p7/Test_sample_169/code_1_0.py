import numpy as np

def model_latitudinal_gradient():
    """
    This script models and explains the hypothesized relationship between latitude 
    and plant Volatile Organic Compound (VOC) diversity based on the 
    'information arms-race' theory.
    """

    # --- Model Parameters ---
    # The effect of latitude is modeled using a linear equation:
    # Diversity = Intercept - (Slope * Latitude)
    # A positive slope value in our equation signifies a negative relationship,
    # as diversity decreases when latitude increases.

    # Alpha (α) Diversity - within a single plant
    # Higher baseline diversity at the equator (latitude=0)
    alpha_intercept = 4.8 
    # The amount diversity decreases for each degree increase in latitude
    alpha_slope_effect = 0.055

    # Beta (β) Diversity - among plants within a site
    # Higher baseline dissimilarity at the equator (latitude=0)
    beta_intercept = 0.9
    # The amount dissimilarity decreases for each degree increase in latitude
    beta_slope_effect = 0.012

    # --- Explanation and Output ---
    print("Ecological Theory: The 'information arms-race' between plants and their parasites is strongest in the tropics (low latitude) and weakens towards the poles (high latitude).")
    print("-" * 70)
    
    # Alpha Diversity Explanation
    print("Prediction for Alpha (α) Diversity (VOCs within a single plant):")
    print("1. In the tropics, intense parasite pressure selects for plants with a more complex chemical arsenal.")
    print("2. Therefore, α-diversity is highest at the equator and decreases as latitude increases.")
    print("3. This represents a NEGATIVE effect of latitude on α-diversity.")
    print("\n   Modeled Equation:")
    print(f"   α-diversity = {alpha_intercept} - {alpha_slope_effect} * (Latitude)")
    print("-" * 70)

    # Beta Diversity Explanation
    print("Prediction for Beta (β) Diversity (VOC differences among plants):")
    print("1. In the tropics, this same pressure selects for high chemical variation among plants to act as a 'moving target' for parasites.")
    print("2. Therefore, β-diversity (dissimilarity) is highest at the equator and decreases as latitude increases.")
    print("3. This also represents a NEGATIVE effect of latitude on β-diversity.")
    print("\n   Modeled Equation:")
    print(f"   β-diversity = {beta_intercept} - {beta_slope_effect} * (Latitude)")
    print("-" * 70)

    print("\nCONCLUSION: Both α-diversity and β-diversity are expected to have a negative relationship with latitude.")

# Execute the function to print the explanation.
model_latitudinal_gradient()
<<<B>>>