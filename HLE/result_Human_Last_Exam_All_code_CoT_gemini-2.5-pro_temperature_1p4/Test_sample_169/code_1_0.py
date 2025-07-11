def model_and_explain_diversity_gradient():
    """
    This script models the theoretical relationship between latitude and plant
    Volatile Organic Compound (VOC) diversity based on ecological principles.

    Principle 1: Biotic pressure (from herbivores, pathogens) is highest at the
                 equator (0 degrees latitude) and decreases towards the poles.
    Principle 2: This intense pressure drives the evolution of both within-plant
                 (alpha) and between-plant (beta) chemical diversity as a defense.

    The model will show that as latitude increases, both diversity metrics are expected
    to decrease, demonstrating a negative effect.
    """
    print("Modeling the effect of Latitude on VOC diversity...")

    # Define a simple linear model where diversity decreases with latitude.
    # The exact numbers are for illustration, but the negative slope is the key concept.

    # Model for Alpha Diversity (within a single plant)
    # Let's assume a maximum Shannon-Weiner index of 4.5 at the equator.
    base_alpha = 4.5
    # And it decreases by 0.06 for every degree of latitude north.
    slope_alpha = -0.06
    equation_alpha = f"Alpha Diversity = {base_alpha} + ({slope_alpha} * Latitude)"

    # Model for Beta Diversity (turnover between plants at a site)
    # Let's assume a maximum beta diversity of 0.9 at the equator.
    base_beta = 0.9
    # And it decreases by 0.012 for every degree of latitude north.
    slope_beta = -0.012
    equation_beta = f"Beta Diversity = {base_beta} + ({slope_beta} * Latitude)"

    print(f"\nModel for Alpha Diversity: {equation_alpha}")
    print(f"Model for Beta Diversity:  {equation_beta}\n")

    # --- Calculations for specific latitudes ---
    latitudes_to_test = [0, 30, 60]

    print("="*60)
    print("Latitude (N)\t| Predicted Alpha Diversity | Predicted Beta Diversity")
    print("="*60)

    for lat in latitudes_to_test:
        # Calculate alpha diversity: alpha = 4.5 + (-0.06 * lat)
        alpha_div = base_alpha + slope_alpha * lat
        # Calculate beta diversity: beta = 0.9 + (-0.012 * lat)
        beta_div = base_beta + slope_beta * lat

        # We use max(0, ...) to ensure diversity doesn't become negative.
        print(f"{lat:<15.1f}\t| {max(0, alpha_div):<25.2f} | {max(0, beta_div):<25.2f}")

    print("="*60)
    print("\nConclusion:")
    print("As latitude increases from the equator (0) to 60N, the models predict that both")
    print("alpha diversity and beta diversity will decrease.")
    print("This indicates a NEGATIVE effect of latitude on both alpha and beta diversity.")

# Run the modeling and explanation function.
model_and_explain_diversity_gradient()