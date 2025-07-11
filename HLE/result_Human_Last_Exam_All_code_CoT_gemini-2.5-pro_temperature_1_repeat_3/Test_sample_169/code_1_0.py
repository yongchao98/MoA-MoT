import numpy as np
from scipy.stats import entropy
from scipy.spatial.distance import braycurtis

def simulate_and_analyze():
    """
    Simulates and analyzes plant VOC diversity at high and low latitude sites
    to illustrate the effect of latitude on alpha and beta diversity.
    """
    # --- 1. Simulation Parameters ---
    high_lat_deg = 60.0
    low_lat_deg = 10.0
    num_plants_per_site = 5
    total_voc_pool_size = 25 # Total possible VOCs in the simulation

    # --- 2. Simulate High Latitude Site (e.g., 60N) ---
    # Expectation: Lower alpha diversity and lower beta diversity.
    high_lat_plants = []
    for _ in range(num_plants_per_site):
        # Lower alpha: each plant expresses fewer VOCs (e.g., 4 to 7)
        num_vocs = np.random.randint(4, 8)
        # Lower beta: VOCs are drawn from a smaller, common pool (e.g., the first 10)
        voc_indices = np.random.choice(range(10), size=num_vocs, replace=False)
        profile = np.zeros(total_voc_pool_size)
        # Assign random relative abundances
        profile[voc_indices] = np.random.uniform(0.1, 1.0, size=num_vocs)
        profile /= profile.sum() # Normalize to get proportions (p_i)
        high_lat_plants.append(profile)

    # --- 3. Simulate Low Latitude Site (e.g., 10N) ---
    # Expectation: Higher alpha diversity and higher beta diversity.
    low_lat_plants = []
    for _ in range(num_plants_per_site):
        # Higher alpha: each plant expresses more VOCs (e.g., 12 to 18)
        num_vocs = np.random.randint(12, 19)
        # Higher beta: VOCs are drawn from the entire, larger pool
        voc_indices = np.random.choice(range(total_voc_pool_size), size=num_vocs, replace=False)
        profile = np.zeros(total_voc_pool_size)
        profile[voc_indices] = np.random.uniform(0.1, 1.0, size=num_vocs)
        profile /= profile.sum() # Normalize
        low_lat_plants.append(profile)

    # --- 4. Calculate and Print Alpha Diversity (Shannon-Wiener Index, H) ---
    # H = -sum(p_i * log(p_i)). We use scipy.stats.entropy which calculates this.
    avg_alpha_high = np.mean([entropy(p) for p in high_lat_plants])
    avg_alpha_low = np.mean([entropy(p) for p in low_lat_plants])

    print("--- Analysis of α-Diversity (Within-Plant Diversity) ---")
    print(f"Average α-diversity at {high_lat_deg}°N: {avg_alpha_high:.4f}")
    print(f"Average α-diversity at {low_lat_deg}°N: {avg_alpha_low:.4f}")
    print("Observation: As latitude increases, α-diversity decreases (Negative relationship).\n")

    # --- 5. Calculate and Print Beta Diversity (Bray-Curtis Dissimilarity) ---
    # We use average pairwise Bray-Curtis dissimilarity as a measure of beta diversity.
    dissimilarities_high = [braycurtis(p1, p2) for i, p1 in enumerate(high_lat_plants) for p2 in high_lat_plants[i+1:]]
    avg_beta_high = np.mean(dissimilarities_high)
    dissimilarities_low = [braycurtis(p1, p2) for i, p1 in enumerate(low_lat_plants) for p2 in low_lat_plants[i+1:]]
    avg_beta_low = np.mean(dissimilarities_low)

    print("--- Analysis of β-Diversity (Among-Plant Variation) ---")
    print(f"Average β-diversity (Bray-Curtis) at {high_lat_deg}°N: {avg_beta_high:.4f}")
    print(f"Average β-diversity (Bray-Curtis) at {low_lat_deg}°N: {avg_beta_low:.4f}")
    print("Observation: As latitude increases, β-diversity decreases (Negative relationship).\n")

    # --- 6. Deriving Illustrative Linear Equations ---
    # Equation form: Diversity = m * Latitude + b
    # Calculate slope (m) and intercept (b) from our two simulated points.
    m_alpha = (avg_alpha_high - avg_alpha_low) / (high_lat_deg - low_lat_deg)
    b_alpha = avg_alpha_high - m_alpha * high_lat_deg

    m_beta = (avg_beta_high - avg_beta_low) / (high_lat_deg - low_lat_deg)
    b_beta = avg_beta_high - m_beta * high_lat_deg

    print("--- Illustrative Linear Equation from Simulation ---")
    print("Based on the two simulated points, the relationship is:")
    print("\nFinal equation for α-Diversity:")
    print(f"α-Diversity = {m_alpha} * Latitude + {b_alpha}")

    print("\nFinal equation for β-Diversity:")
    print(f"β-Diversity = {m_beta} * Latitude + {b_beta}")

# Run the simulation and analysis
simulate_and_analyze()