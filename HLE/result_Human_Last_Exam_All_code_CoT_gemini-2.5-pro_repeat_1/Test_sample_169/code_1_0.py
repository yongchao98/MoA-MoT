import pandas as pd
import numpy as np
import statsmodels.api as sm

def analyze_voc_diversity_gradient():
    """
    Analyzes the theoretical effect of latitude on plant VOC diversity
    based on the biotic arms-race hypothesis.
    """
    # --- Step 1: Explanation of the Ecological Hypothesis ---
    # The Latitudinal Diversity Gradient (LDG) shows that biodiversity is highest in the tropics
    # (low latitudes) and decreases towards the poles. A key driver is the "biotic arms-race"
    # between plants and their parasites (herbivores, pathogens), which is more intense in the tropics.
    #
    # - Alpha (α) diversity (VOCs within a single plant): In the tropics, higher parasite pressure
    #   forces individual plants to evolve a more complex and diverse suite of VOCs for defense.
    #   Therefore, α-diversity increases as latitude decreases.
    #
    # - Beta (β) diversity (VOC variation among plants at a site): To prevent parasites from
    #   easily adapting and spreading, there is strong selection for chemical differentiation among
    #   plants in the tropics. This high variation in VOC profiles among individuals leads to high
    #   β-diversity. Therefore, β-diversity also increases as latitude decreases.
    #
    # Conclusion: The effect of latitude on both α and β diversity is negative.
    # As the latitude value increases (0 -> 60), diversity decreases.

    # --- Step 2: Simulate Data to Demonstrate the Hypothesis ---
    # We create a synthetic dataset reflecting these expected negative relationships.
    np.random.seed(42)  # for reproducible results
    num_sites = 100
    
    # Latitudes are sampled uniformly from the equator (0) to 60 degrees.
    latitude = np.linspace(0, 60, num_sites)

    # Generate alpha diversity with a clear negative trend plus some random noise.
    # We assume a baseline diversity of 3.5 at the equator, which decreases with latitude.
    # Equation: alpha_diversity = 3.5 - 0.03 * latitude + noise
    alpha_diversity = 3.5 - 0.03 * latitude + np.random.normal(loc=0, scale=0.2, size=num_sites)

    # Generate beta diversity with a clear negative trend plus some random noise.
    # We assume a baseline diversity of 0.8 at the equator, which decreases with latitude.
    # Equation: beta_diversity = 0.8 - 0.01 * latitude + noise
    beta_diversity = 0.8 - 0.01 * latitude + np.random.normal(loc=0, scale=0.05, size=num_sites)

    # --- Step 3: Perform Linear Regression Analysis ---
    # We will use ordinary least squares (OLS) regression to find the equation
    # describing the relationship between latitude and diversity.
    
    # Add a constant (for the intercept) to our predictor variable
    X = sm.add_constant(latitude)

    # Fit the model for Alpha Diversity
    model_alpha = sm.OLS(alpha_diversity, X).fit()
    intercept_alpha, slope_alpha = model_alpha.params

    # Fit the model for Beta Diversity
    model_beta = sm.OLS(beta_diversity, X).fit()
    intercept_beta, slope_beta = model_beta.params

    # --- Step 4: Output the Results and Final Equations ---
    print("This script simulates and analyzes the effect of latitude on plant VOC diversity.")
    print("Based on the biotic arms-race hypothesis, we expect a negative relationship for both alpha and beta diversity.\n")

    print("--- Analysis of Alpha (α) Diversity vs. Latitude ---")
    print(f"The calculated effect (slope) of latitude on α-diversity is {slope_alpha:.4f}.")
    print("The final equation is:")
    print(f"α-diversity = {intercept_alpha:.4f} + ({slope_alpha:.4f} * latitude)")
    print("\n")

    print("--- Analysis of Beta (β) Diversity vs. Latitude ---")
    print(f"The calculated effect (slope) of latitude on β-diversity is {slope_beta:.4f}.")
    print("The final equation is:")
    print(f"β-diversity = {intercept_beta:.4f} + ({slope_beta:.4f} * latitude)")
    print("\n")

    print("Conclusion: Since both slopes are negative, the analysis supports a negative effect of latitude on both α and β diversity.")

if __name__ == '__main__':
    analyze_voc_diversity_gradient()
<<<B>>>