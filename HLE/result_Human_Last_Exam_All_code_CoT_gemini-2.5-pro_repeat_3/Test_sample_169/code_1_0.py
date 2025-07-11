def model_voc_diversity_gradient():
    """
    Models and explains the effect of latitude on plant VOC diversity based on biotic pressures.
    """
    # --- Model Parameters ---
    # These parameters define the conceptual linear relationship: Diversity = Intercept + Slope * Latitude
    # We choose an intercept (diversity at the equator, latitude=0) and a negative slope.

    # Alpha Diversity (within-plant)
    alpha_intercept = 4.5  # Assumed Shannon-Weiner index at the equator
    alpha_slope = -0.05    # Negative slope: diversity decreases as latitude increases

    # Beta Diversity (between-plant, within-site)
    beta_intercept = 0.7   # Assumed beta diversity value at the equator
    beta_slope = -0.01     # Negative slope: diversity decreases as latitude increases

    # --- Latitudes of Interest ---
    low_latitude = 0   # Equator
    high_latitude = 60 # Northern latitude

    # --- Calculate Diversity at each Latitude ---
    predicted_alpha_low_lat = alpha_intercept + alpha_slope * low_latitude
    predicted_alpha_high_lat = alpha_intercept + alpha_slope * high_latitude

    predicted_beta_low_lat = beta_intercept + beta_slope * low_latitude
    predicted_beta_high_lat = beta_intercept + beta_slope * high_latitude

    # --- Print Explanation and Results ---
    print("Ecological Rationale:")
    print("Biotic pressure from parasites is highest at low latitudes (tropics) and decreases towards the poles.")
    print("This intense pressure drives the evolution of greater chemical diversity in plants.\n")

    print("--- Alpha (α) Diversity: Within a single plant ---")
    print("Prediction: Higher parasite pressure at low latitudes selects for more complex chemical defenses within each plant.")
    print(f"Model Equation: α-diversity = {alpha_intercept} + ({alpha_slope} * Latitude)")
    print(f"Predicted α-diversity at {low_latitude}° (equator): {predicted_alpha_low_lat:.2f}")
    print(f"Predicted α-diversity at {high_latitude}°: {predicted_alpha_high_lat:.2f}")
    print("Result: The relationship between latitude and α-diversity is NEGATIVE.\n")

    print("--- Beta (β) Diversity: Variation among plants at a site ---")
    print("Prediction: High pressure from specialized parasites at low latitudes favors chemical uniqueness among neighboring plants, increasing community-wide chemical variation.")
    print(f"Model Equation: β-diversity = {beta_intercept} + ({beta_slope} * Latitude)")
    print(f"Predicted β-diversity at {low_latitude}° (equator): {predicted_beta_low_lat:.2f}")
    print(f"Predicted β-diversity at {high_latitude}°: {predicted_beta_high_lat:.2f}")
    print("Result: The relationship between latitude and β-diversity is NEGATIVE.")

# Run the model
model_voc_diversity_gradient()