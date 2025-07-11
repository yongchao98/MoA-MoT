def demonstrate_latitudinal_gradient(latitude):
    """
    This function provides a conceptual model of how latitude affects
    plant chemical diversity based on the intensity of biotic interactions.

    Args:
        latitude (int): The latitude in degrees (0 = equator, 90 = pole).
    """

    # --- Model Parameters ---
    # Baseline diversity at the equator (latitude = 0), representing the most intense biotic interactions.
    # These are arbitrary values for demonstration.
    base_alpha_diversity = 3.5  # e.g., Shannon index for VOCs within a plant
    base_beta_diversity = 0.8   # e.g., Average Bray-Curtis dissimilarity among plants

    # Effect coefficients. A negative value means diversity decreases as latitude increases.
    # This simulates the weakening of the plant-parasite arms race away from the tropics.
    alpha_effect_per_degree = -0.04
    beta_effect_per_degree = -0.01

    # --- Calculation ---
    # The equation models diversity as a function of a baseline at the equator
    # plus the cumulative negative effect of increasing latitude.
    final_alpha = base_alpha_diversity + alpha_effect_per_degree * latitude
    final_beta = base_beta_diversity + beta_effect_per_degree * latitude
    
    # Ensure diversity doesn't fall below a realistic minimum
    final_alpha = max(0.5, final_alpha)
    final_beta = max(0.1, final_beta)

    print(f"--- Calculating Diversity for Latitude: {latitude}°N ---")
    print("Alpha-diversity (within a single plant):")
    print(f"Equation: {base_alpha_diversity} + ({alpha_effect_per_degree} * {latitude}) = {final_alpha:.2f}")
    print("Beta-diversity (dissimilarity among plants):")
    print(f"Equation: {base_beta_diversity} + ({beta_effect_per_degree} * {latitude}) = {final_beta:.2f}\n")


# --- Explanation ---
print("Hypothesis: The plant-parasite 'arms race' is most intense in the tropics (low latitude) and weakens towards the poles (high latitude).\n")
print("1. Alpha-diversity: Higher parasite pressure in the tropics forces individual plants to evolve a more diverse chemical arsenal. Thus, as latitude INCREASES, alpha-diversity DECREASES (a negative effect).")
print("2. Beta-diversity: Higher parasite pressure and species richness in the tropics favor chemical uniqueness among neighboring plants to avoid specialized enemies. Thus, as latitude INCREASES, the dissimilarity among plants DECREASES (a negative effect).\n")

# --- Demonstration ---
print("--- Model Demonstration ---")
print("Comparing a low-latitude (tropical) site with a high-latitude (temperate) site:\n")
# A low-latitude site (closer to the equator)
demonstrate_latitudinal_gradient(latitude=10)

# A high-latitude site (further from the equator)
demonstrate_latitudinal_gradient(latitude=60)
