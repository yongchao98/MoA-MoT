import math

def model_diversity_gradient():
    """
    This script models the relationship between latitude and plant VOC diversity
    based on the principle of the latitudinal gradient in biotic interaction strength.
    """
    print("Modeling the effect of latitude on plant VOC alpha and beta diversity.")
    print("Assumption: Biotic pressure (from parasites/herbivores) increases as latitude decreases.")
    print("-" * 70)
    print(f"{'Latitude (N)':<15} | {'Relative Pressure':<20} | {'Alpha Diversity (Simulated)':<28} | {'Beta Diversity (Simulated)':<25}")
    print("-" * 70)

    # We will model latitudes from 60 degrees North down to the Equator (0 degrees).
    # We use a step of -10 to show the trend clearly.
    latitudes = range(60, -1, -10)

    # Constants for the simulation. These are arbitrary values to show the relationship.
    # We use a non-linear relationship to model that pressure increases sharply in the tropics.
    PRESSURE_CONSTANT = 500.0
    ALPHA_DIVERSITY_FACTOR = 0.5
    BETA_DIVERSITY_FACTOR = 0.8

    for lat in latitudes:
        # Model biotic pressure as inversely related to latitude.
        # Adding 1 to latitude avoids division by zero at the equator.
        # The exponent makes the increase in pressure more pronounced at low latitudes.
        relative_pressure = PRESSURE_CONSTANT / ((lat + 1)**0.5)

        # Model alpha and beta diversity as being positively correlated with biotic pressure.
        # A plant needs more compounds (alpha) to fight more enemies.
        # A population needs more variation (beta) to prevent epidemics.
        alpha_diversity = ALPHA_DIVERSITY_FACTOR * relative_pressure
        beta_diversity = BETA_DIVERSITY_FACTOR * relative_pressure

        # The final equation for alpha diversity for a given latitude 'lat' is:
        # alpha_diversity = ALPHA_DIVERSITY_FACTOR * (PRESSURE_CONSTANT / ((lat + 1)**0.5))
        # Example for Latitude = 60: alpha_diversity = 0.5 * (500.0 / ((60 + 1)**0.5))
        # We will print the final values for each step.
        print(f"{lat:<15} | {relative_pressure:<20.2f} | {alpha_diversity:<28.2f} | {beta_diversity:<25.2f}")

    print("-" * 70)
    print("\nConclusion from model:")
    print("As Latitude decreases (moving towards the equator), both Alpha and Beta diversity increase.")
    print("This demonstrates a 'negative' effect, as diversity decreases when latitude increases.")

# Execute the simulation
model_diversity_gradient()
