import pandas as pd

def model_the_relationship():
    """
    This function models the theoretical relationship between latitude and
    plant VOC diversity based on the principle of the latitudinal biotic pressure gradient.
    """
    print("Plan: This model will demonstrate the expected effect of latitude on plant VOC diversity.")
    print("1. The model assumes that biotic pressure (from parasites, etc.) is the primary driver.")
    print("2. Biotic pressure is highest at the equator (0 degrees) and decreases towards the poles.")
    print("3. Both within-plant diversity (alpha) and between-plant diversity (beta) are positively correlated with biotic pressure.")
    print("4. Therefore, both alpha and beta diversity should show a negative relationship with latitude.\n")

    # We will model a few points from the study range (equator to 60 degrees N)
    latitudes = [0, 20, 40, 60]

    # Model parameters (arbitrary values to show the relationship)
    # Let's say max alpha diversity (at equator) is 4.5 and max beta diversity is 0.8
    max_alpha = 4.5
    max_beta = 0.8
    # We model diversity as decreasing linearly up to the pole (90 degrees)
    pole_latitude = 90.0

    results = []
    for lat in latitudes:
        # Calculate the diversity as a function of latitude.
        # The term (1 - lat / pole_latitude) decreases from 1 to 0 as latitude increases.
        # This equation models the negative relationship.
        alpha_diversity = max_alpha * (1 - lat / pole_latitude)
        beta_diversity = max_beta * (1 - lat / pole_latitude)
        results.append({
            "Latitude (Degrees)": lat,
            "Modeled Alpha Diversity": f"{alpha_diversity:.2f}",
            "Modeled Beta Diversity": f"{beta_diversity:.2f}"
        })
        
    # Using pandas for a clean, tabular output
    df = pd.DataFrame(results)
    print("--- Illustrative Model Output ---")
    print(df.to_string(index=False))
    print("\n--- Conclusion ---")
    print("As shown in the model, when Latitude increases, both Alpha and Beta diversity decrease.")
    print("This indicates a 'negative' direction of effect for both.")

model_the_relationship()