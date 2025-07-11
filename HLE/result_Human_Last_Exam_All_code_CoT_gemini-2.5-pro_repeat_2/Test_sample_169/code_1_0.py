import numpy as np
from scipy import stats

def demonstrate_voc_latitude_gradient():
    """
    This function simulates and demonstrates the theoretical relationship
    between latitude and plant VOC diversity (alpha and beta).
    """
    # --- 1. Setup Simulation Parameters ---
    # We simulate sites from the equator (0) to 60 degrees North.
    np.random.seed(42) # for reproducible results
    latitudes = np.linspace(0, 60, 30) # 30 sites from 0 to 60 degrees

    # --- 2. Model the Ecological Principle ---
    # Assumption: Biotic pressure (from herbivores/parasites) is highest at the
    # equator (latitude 0) and decreases linearly towards the poles.
    # We model it as a value from 10 (max pressure) down to some lower value.
    max_pressure = 10
    biotic_pressure = max_pressure - (latitudes / 60) * 8

    # --- 3. Generate Simulated Diversity Data ---
    # Assumption: Both alpha and beta diversity are positively correlated with
    # biotic pressure, but with some natural random variation (noise).
    # alpha_diversity: Complexity of VOCs in a single plant.
    # beta_diversity: Variation of VOCs among plants at a site.
    alpha_noise = np.random.normal(0, 0.8, size=latitudes.shape)
    beta_noise = np.random.normal(0, 1.0, size=latitudes.shape)

    simulated_alpha_diversity = biotic_pressure + alpha_noise
    simulated_beta_diversity = biotic_pressure * 0.7 + beta_noise # Beta diversity might have a different scale

    # Ensure diversity values are not negative
    simulated_alpha_diversity[simulated_alpha_diversity < 0] = 0
    simulated_beta_diversity[simulated_beta_diversity < 0] = 0

    # --- 4. Analyze the Relationship using Linear Regression ---
    # We check if the slope of the relationship between latitude and diversity is negative.
    alpha_reg = stats.linregress(latitudes, simulated_alpha_diversity)
    beta_reg = stats.linregress(latitudes, simulated_beta_diversity)

    # --- 5. Print the Results and Explanation ---
    print("--- Ecological Theory Simulation: Latitude vs. Plant VOC Diversity ---")
    print("\nBased on the principle that biotic pressure decreases from the equator to the poles,")
    print("we expect both the complexity (alpha-diversity) and site-level variation (beta-diversity)")
    print("of plant chemical defenses to decrease with increasing latitude.")
    print("-" * 70)

    # Alpha Diversity Results
    print("\n1. Analysis of α-diversity (VOC complexity within a single plant):")
    print(f"The linear regression equation is: alpha_diversity = {alpha_reg.slope:.4f} * latitude + {alpha_reg.intercept:.4f}")
    print(f"The slope is {alpha_reg.slope:.4f}. Since this value is negative, it indicates that as latitude increases,")
    print("the predicted alpha-diversity of VOCs decreases.")
    print("-" * 70)

    # Beta Diversity Results
    print("\n2. Analysis of β-diversity (VOC variation among plants at a site):")
    print(f"The linear regression equation is: beta_diversity = {beta_reg.slope:.4f} * latitude + {beta_reg.intercept:.4f}")
    print(f"The slope is {beta_reg.slope:.4f}. Since this value is also negative, it indicates that as latitude increases,")
    print("the predicted beta-diversity of VOCs decreases.")
    print("-" * 70)

    print("\nConclusion: The simulation supports the theoretical expectation of a negative relationship for both alpha and beta diversity.")

if __name__ == '__main__':
    demonstrate_voc_latitude_gradient()
<<<B>>>