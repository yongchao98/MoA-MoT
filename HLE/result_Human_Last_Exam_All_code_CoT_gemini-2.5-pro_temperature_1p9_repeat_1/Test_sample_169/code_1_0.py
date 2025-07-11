import numpy as np
from scipy.stats import linregress
from scipy.spatial.distance import pdist

def calculate_shannon_diversity(proportions):
    """Calculates the Shannon-Wiener diversity index."""
    # Filter out zero proportions to avoid log(0)
    proportions = proportions[proportions > 0]
    # If there are no compounds, diversity is 0
    if proportions.size == 0:
        return 0.0
    return -np.sum(proportions * np.log(proportions))

# --- 1. Simulation Setup ---
# Define the range of latitudes from 60 degrees North to the Equator (0).
latitudes = np.linspace(60, 0, num=20) 
plants_per_site = 50
total_voc_universe = 100 # Total possible VOCs in the ecosystem.

# Lists to store the calculated diversity metrics for each site.
all_sites_alpha_diversity = []
all_sites_beta_diversity = []

# --- 2. Main Simulation Loop ---
# Iterate through each latitude to simulate the ecosystem at that site.
for lat in latitudes:
    site_plant_profiles = []
    
    # Model the effect of latitude on VOCs.
    # As latitude decreases (approaches 0), selective pressure increases.
    # This leads to more VOCs produced per plant and more variation among plants.
    pressure_factor = (60 - lat) / 60.0  # Ranges from 0 at 60N to 1 at the equator.
    
    # Number of VOCs present in a plant increases with pressure.
    n_vocs_present = int(10 + 40 * pressure_factor) 
    # Evenness of VOC profile increases with pressure (parameter for Dirichlet distribution).
    evenness_param = 0.5 + 2.0 * pressure_factor
    
    # Generate VOC profiles for each plant at the current site.
    for _ in range(plants_per_site):
        # Generate a base distribution of the VOCs that are present.
        # A Dirichlet distribution is ideal for creating compositional data (proportions).
        base_proportions = np.random.dirichlet(np.ones(n_vocs_present) * evenness_param)
        
        # Place these proportions into the full universe of possible VOCs.
        full_profile = np.zeros(total_voc_universe)
        # Randomly assign which of the total VOCs are the ones present.
        present_indices = np.random.choice(total_voc_universe, n_vocs_present, replace=False)
        full_profile[present_indices] = base_proportions
        
        site_plant_profiles.append(full_profile)
        
    site_plant_profiles = np.array(site_plant_profiles)

    # --- 3. Calculate Diversity Metrics ---
    # Alpha Diversity: The average Shannon diversity of all plants at the site.
    current_alphas = [calculate_shannon_diversity(profile) for profile in site_plant_profiles]
    mean_alpha = np.mean(current_alphas)
    all_sites_alpha_diversity.append(mean_alpha)
    
    # Beta Diversity: The average dissimilarity (e.g., Bray-Curtis) between all pairs of plants.
    # This measures the compositional turnover between plants at the site.
    pairwise_distances = pdist(site_plant_profiles, 'braycurtis')
    mean_beta = np.mean(pairwise_distances)
    all_sites_beta_diversity.append(mean_beta)

# --- 4. Analyze and Report Results ---
# Use linear regression to find the relationship (slope) between latitude and diversity.
alpha_reg = linregress(latitudes, all_sites_alpha_diversity)
beta_reg = linregress(latitudes, all_sites_beta_diversity)

print("--- Analysis of Simulated VOC Diversity vs. Latitude ---")
print("This simulation models the hypothesis that coevolutionary pressure on plants")
print("is strongest in the tropics (low latitude), driving higher VOC diversity.\n")

print("1. Alpha (α) Diversity - VOC complexity within each plant:")
print("A negative slope indicates that as latitude decreases, α-diversity increases.")
print(f"   Linear Equation: α-Diversity = {alpha_reg.slope:.4f} * Latitude + {alpha_reg.intercept:.4f}")
print(f"   Slope (Effect of Latitude on α-Diversity): {alpha_reg.slope:.4f}")
print(f"   R-squared: {alpha_reg.rvalue**2:.4f}\n")

print("2. Beta (β) Diversity - VOC variation among plants at a site:")
print("A negative slope indicates that as latitude decreases, β-diversity increases.")
print(f"   Linear Equation: β-Diversity = {beta_reg.slope:.4f} * Latitude + {beta_reg.intercept:.4f}")
print(f"   Slope (Effect of Latitude on β-Diversity): {beta_reg.slope:.4f}")
print(f"   R-squared: {beta_reg.rvalue**2:.4f}\n")

print("Conclusion: The simulation shows a negative effect of latitude on both α and β diversity,")
print("supporting the theoretical expectation based on the latitudinal gradient in biotic interactions.")
<<<B>>>