import sys

def model_voc_diversity_gradient():
    """
    Models the hypothetical relationship between latitude and plant VOC diversity
    based on the biotic arms-race hypothesis.
    """
    # --- Model Parameters ---
    # These parameters are chosen to illustrate the predicted trend.

    # Alpha Diversity (within a single plant)
    # Assumes high diversity at the equator (lat=0) that decreases with latitude.
    alpha_intercept = 3.8  # Max Shannon-Weiner diversity at the equator
    alpha_slope = 0.04     # Rate of diversity loss per degree of latitude

    # Beta Diversity (variation among plants at a site)
    # Assumes high variation at the equator that decreases with latitude.
    beta_intercept = 0.9   # Max beta-diversity (e.g., Bray-Curtis) at the equator
    beta_slope = 0.012     # Rate of variation loss per degree of latitude

    # Latitudes to examine (Equator, Subtropics, High Latitude)
    latitudes = [0, 30, 60]

    # --- Output ---
    print("--- Modeling the Effect of Latitude on Plant VOC Diversity ---")
    print("\nHypothesis: Both alpha and beta diversity decrease with increasing latitude due to reduced enemy pressure.")

    print(f"\nAlpha Diversity Model: Diversity = {alpha_intercept} - ({alpha_slope} * Latitude)")
    print(f"Beta Diversity Model:  Diversity = {beta_intercept} - ({beta_slope} * Latitude)")
    print("-" * 60)

    for lat in latitudes:
        # Calculate predicted alpha diversity
        predicted_alpha = alpha_intercept - (alpha_slope * lat)
        
        # Calculate predicted beta diversity
        predicted_beta = beta_intercept - (beta_slope * lat)

        print(f"\nCalculations for Latitude = {lat}° N:")
        
        # Print the full equation and result for alpha diversity
        sys.stdout.write(f"  Predicted Alpha Diversity = {alpha_intercept} - ({alpha_slope} * {lat}) = ")
        sys.stdout.flush()
        print(f"{predicted_alpha:.3f}")

        # Print the full equation and result for beta diversity
        sys.stdout.write(f"  Predicted Beta Diversity  = {beta_intercept} - ({beta_slope} * {lat}) = ")
        sys.stdout.flush()
        print(f"{predicted_beta:.3f}")

    print("-" * 60)
    print("Conclusion: The model shows a negative effect of latitude on both alpha and beta diversity.")

if __name__ == '__main__':
    model_voc_diversity_gradient()