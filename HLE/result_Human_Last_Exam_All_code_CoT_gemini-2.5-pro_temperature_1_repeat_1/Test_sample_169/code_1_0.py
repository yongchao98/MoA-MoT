import numpy as np

def simulate_voc_diversity():
    """
    Simulates and prints hypothetical plant VOC diversity data along a latitudinal gradient.
    This code demonstrates the predicted negative relationship between latitude and
    both alpha (within-plant) and beta (among-plant) diversity of VOCs.
    """
    print("This simulation demonstrates the expected ecological patterns based on the 'arms-race' hypothesis.")
    print("As biotic pressure from parasites is higher in the tropics (low latitude), both within-plant (alpha) and among-plant (beta) VOC diversity are expected to be higher.\n")

    # Define the parameters for our linear model.
    # The negative slope signifies a negative relationship with latitude.
    # Equation form: Diversity = Intercept - Slope * Latitude + Noise

    # Parameters for alpha diversity (e.g., Shannon-Weiner Index)
    alpha_intercept = 4.5
    alpha_slope = 0.05
    
    # Parameters for beta diversity (e.g., a dissimilarity index from 0 to 1)
    beta_intercept = 0.9
    beta_slope = 0.012

    # Print the model equations
    print("--- Model Equations ---")
    print(f"Alpha Diversity = {alpha_intercept} - {alpha_slope} * Latitude")
    print(f"Beta Diversity  = {beta_intercept} - {beta_slope} * Latitude")
    print("-----------------------\n")
    
    # Latitudes from the study design (equator to 60 degrees north)
    latitudes_to_sample = [0, 15, 30, 45, 60]

    print("--- Simulated Data ---")
    print(f"{'Latitude (°N)':<15}{'Predicted α-Diversity':<25}{'Predicted β-Diversity':<25}")
    print("-" * 65)

    # Set a seed for reproducibility of the random noise
    np.random.seed(42)

    for lat in latitudes_to_sample:
        # Add a small amount of random noise to make it look more realistic
        noise_alpha = np.random.normal(0, 0.1)
        noise_beta = np.random.normal(0, 0.02)
        
        # Calculate the diversity values using the linear model
        # We use max(..., 0) to ensure diversity isn't negative
        alpha_diversity = max(alpha_intercept - alpha_slope * lat + noise_alpha, 0)
        beta_diversity = max(beta_intercept - beta_slope * lat + noise_beta, 0)
        
        print(f"{lat:<15.1f}{alpha_diversity:<25.3f}{beta_diversity:<25.3f}")

# Run the simulation
simulate_voc_diversity()