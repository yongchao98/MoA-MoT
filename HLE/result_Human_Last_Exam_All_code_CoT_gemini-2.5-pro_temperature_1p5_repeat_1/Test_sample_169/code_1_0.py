def solve_voc_diversity_problem():
    """
    This script provides a conceptual model to determine the effect of latitude on
    plant VOC diversity based on the "information arms-race" hypothesis.
    """

    # Step 1: Model the relationship between latitude and biotic pressure.
    # The "information arms-race" (biotic pressure) is most intense at low latitudes (tropics)
    # and decreases towards higher latitudes.
    def get_biotic_pressure(latitude):
        """A simple model where pressure is inversely related to latitude."""
        # We use a simple linear model for illustration. Max pressure at lat 0.
        return max(0, 100 - 1.5 * latitude)

    # Step 2: Model how biotic pressure influences alpha and beta diversity.
    # Higher pressure demands more complex defenses, increasing both within-plant (alpha)
    # and between-plant (beta) diversity.
    def get_alpha_diversity(pressure):
        """Alpha diversity (within a plant) increases with pressure."""
        return pressure * 1.2

    def get_beta_diversity(pressure):
        """Beta diversity (among plants) increases with pressure."""
        return pressure * 0.8

    # Step 3: Compare a low-latitude site with a high-latitude site.
    low_latitude_site = 0  # Equator
    high_latitude_site = 60 # Temperate/Boreal

    # Calculate values for the low-latitude site
    pressure_low_lat = get_biotic_pressure(low_latitude_site)
    alpha_low_lat = get_alpha_diversity(pressure_low_lat)
    beta_low_lat = get_beta_diversity(pressure_low_lat)

    # Calculate values for the high-latitude site
    pressure_high_lat = get_biotic_pressure(high_latitude_site)
    alpha_high_lat = get_alpha_diversity(pressure_high_lat)
    beta_high_lat = get_beta_diversity(pressure_high_lat)

    # Step 4: Determine the direction of the effect.
    # An effect is "negative" if the value decreases as latitude increases.
    alpha_effect = "negative" if alpha_high_lat < alpha_low_lat else "positive"
    beta_effect = "negative" if beta_high_lat < beta_low_lat else "positive"

    # Step 5: Print the detailed reasoning and conclusion.
    print("Ecological Reasoning Summary:")
    print("-" * 30)
    print("1. Biotic pressure (the 'arms-race') is highest at the equator (low latitude) and decreases towards the poles (high latitude).")
    print(f"   - Modelled Pressure at {low_latitude_site}°N: {pressure_low_lat:.1f}")
    print(f"   - Modelled Pressure at {high_latitude_site}°N: {pressure_high_lat:.1f}\n")

    print("2. Alpha (α) Diversity (VOCs within a single plant):")
    print("   - High pressure selects for more complex individual chemical defenses.")
    print(f"   - Modelled α-Diversity at {low_latitude_site}°N: {alpha_low_lat:.1f}")
    print(f"   - Modelled α-Diversity at {high_latitude_site}°N: {alpha_high_lat:.1f}")
    print(f"   - Result: As latitude increases, α-diversity decreases. The effect is '{alpha_effect}'.\n")

    print("3. Beta (β) Diversity (VOC variation among plants at a site):")
    print("   - High pressure selects for high variation among individuals to prevent parasite specialization.")
    print(f"   - Modelled β-Diversity at {low_latitude_site}°N: {beta_low_lat:.1f}")
    print(f"   - Modelled β-Diversity at {high_latitude_site}°N: {beta_high_lat:.1f}")
    print(f"   - Result: As latitude increases, β-diversity also decreases. The effect is '{beta_effect}'.\n")

    print("Final Conclusion:")
    print(f"The predicted direction of effect is '{alpha_effect}' for α-diversity and '{beta_effect}' for β-diversity.")

# Execute the function to print the analysis.
solve_voc_diversity_problem()