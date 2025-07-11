def model_diversity_gradient():
    """
    Models and explains the effect of latitude on plant VOC diversity
    based on the coevolutionary arms-race hypothesis.
    """
    # --- Step 1: Define the core ecological relationship ---
    # We model pest pressure as being highest at the equator (latitude=0) and
    # decreasing linearly as latitude increases.
    # Equation: Pest_Pressure = p_slope * latitude + p_intercept
    p_intercept = 100  # Relative pest pressure at the equator (latitude 0)
    p_slope = -1.5     # Decrease in pest pressure per degree of latitude increase

    # --- Step 2: Define how diversity scales with pressure ---
    # We model both alpha and beta diversity as being positively correlated
    # with the intensity of pest pressure.
    alpha_scaling_factor = 0.5  # Arbitrary factor for alpha diversity
    beta_scaling_factor = 0.8   # Arbitrary factor for beta diversity

    # --- Step 3: Derive the final equations for diversity vs. latitude ---
    # The relationship is: Diversity = scaling_factor * Pest_Pressure
    # By substituting the pest pressure equation, we get:
    # Diversity = scaling_factor * (p_slope * latitude + p_intercept)
    # This simplifies to: Diversity = (scaling_factor * p_slope) * latitude + (scaling_factor * p_intercept)
    # This is a linear equation of the form: y = m*x + c

    # For alpha diversity:
    alpha_slope = alpha_scaling_factor * p_slope
    alpha_intercept = alpha_scaling_factor * p_intercept

    # For beta diversity:
    beta_slope = beta_scaling_factor * p_slope
    beta_intercept = beta_scaling_factor * p_intercept

    # --- Step 4: Output the explanation and equations ---
    print("Based on the 'arms-race' hypothesis, both alpha and beta VOC diversity are expected to decrease as latitude increases.")
    print("This implies a 'negative' effect of latitude on both diversity metrics.")
    print("\n--- Model Equations (Diversity = slope * Latitude + intercept) ---\n")

    print("Alpha (within-plant) Diversity Equation:")
    print(f"Alpha_Diversity = {alpha_slope:.2f} * Latitude + {alpha_intercept:.2f}")
    print(f"The slope is {alpha_slope:.2f}, which is negative.\n")

    print("Beta (between-plant) Diversity Equation:")
    print(f"Beta_Diversity = {beta_slope:.2f} * Latitude + {beta_intercept:.2f}")
    print(f"The slope is {beta_slope:.2f}, which is also negative.\n")

    # --- Step 5: Show example values ---
    print("--- Example Values ---")
    print("Latitude (deg N) | Alpha Diversity | Beta Diversity")
    print("-------------------------------------------------")
    for lat in [0, 20, 40, 60]:
        alpha_div = alpha_slope * lat + alpha_intercept
        beta_div = beta_slope * lat + beta_intercept
        print(f"{lat:16} | {alpha_div:15.1f} | {beta_div:14.1f}")

model_diversity_gradient()