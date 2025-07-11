import numpy as np
from scipy import stats

def analyze_latitude_diversity_effect():
    """
    Analyzes and simulates the effect of latitude on plant VOC diversity
    based on ecological principles.

    The analysis is based on two key ideas related to the "information arms-race":
    1.  Increased Biotic Pressure at Low Latitudes: The tropics (low latitude) have
        stronger and more consistent pressure from herbivores and parasites than
        temperate zones (high latitude).
    2.  Chemical Defense Response: Plants evolve VOCs as a primary defense.

    α-diversity (within an individual plant): To combat a wider array of enemies in the tropics,
    individual plants are selected to produce a more diverse suite of VOCs.
    Therefore, as latitude decreases, α-diversity increases (a negative relationship).

    β-diversity (among plants at a site): High parasite pressure in the tropics favors
    chemical uniqueness among plants. If individuals in a population have different VOC
    profiles, it is harder for a parasite to adapt and cause an epidemic. This increases
    the chemical turnover, or β-diversity, at the site.
    Therefore, as latitude decreases, β-diversity also increases (a negative relationship).
    """

    # --- Simulation Setup ---
    # We will simulate the relationship to make the concept clear.
    # Latitudes from 60 degrees North down to the Equator (0)
    latitudes = np.linspace(60, 0, 100)

    # --- Model the Relationship ---
    # Model 1: Biotic pressure increases as latitude decreases. We add some random noise.
    # We use a simple linear model: Pressure = Constant - (slope * latitude)
    biotic_pressure = 100 - 1.2 * latitudes + np.random.normal(0, 5, size=100)

    # Model 2: Alpha diversity increases with biotic pressure.
    alpha_diversity = 2.5 + 0.05 * biotic_pressure + np.random.normal(0, 0.2, size=100)

    # Model 3: Beta diversity also increases with biotic pressure.
    beta_diversity = 0.2 + 0.01 * biotic_pressure + np.random.normal(0, 0.05, size=100)

    # --- Analyze the Direction of Effect ---
    # We use linear regression to find the slope between latitude and diversity.
    # A negative slope means diversity decreases as latitude increases.
    alpha_slope, alpha_intercept, _, _, _ = stats.linregress(latitudes, alpha_diversity)
    beta_slope, beta_intercept, _, _, _ = stats.linregress(latitudes, beta_diversity)

    # --- Print the Conclusion ---
    print("Based on ecological principles and the simulation:")
    print("1. Relationship between Latitude and Alpha (α) Diversity:")
    print(f"   The simulated equation is: α-diversity ≈ {alpha_intercept:.3f} + ({alpha_slope:.3f}) * Latitude")
    print("   The slope is negative, indicating α-diversity increases as latitude decreases.")
    print("\n2. Relationship between Latitude and Beta (β) Diversity:")
    print(f"   The simulated equation is: β-diversity ≈ {beta_intercept:.3f} + ({beta_slope:.3f}) * Latitude")
    print("   The slope is also negative, indicating β-diversity increases as latitude decreases.")
    print("\nConclusion: The direction of effect of latitude is negative for both α-diversity and β-diversity.")

# Execute the analysis
analyze_latitude_diversity_effect()