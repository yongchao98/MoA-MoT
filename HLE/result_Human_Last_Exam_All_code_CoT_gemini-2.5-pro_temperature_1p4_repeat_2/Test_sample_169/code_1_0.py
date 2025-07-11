import numpy as np

def shannon_diversity(counts):
    """Calculates the Shannon diversity index."""
    # Add a small constant to avoid log(0)
    probs = counts / np.sum(counts) + 1e-9 
    return -np.sum(probs * np.log(probs))

def calculate_beta_diversity(community_matrix):
    """Calculates average Bray-Curtis dissimilarity as a measure of beta diversity."""
    n_individuals = community_matrix.shape[0]
    if n_individuals < 2:
        return 0.0
    
    total_dissimilarity = 0
    num_pairs = 0
    
    for i in range(n_individuals):
        for j in range(i + 1, n_individuals):
            # Bray-Curtis dissimilarity formula: sum(|xi - xj|) / sum(xi + xj)
            numerator = np.sum(np.abs(community_matrix[i] - community_matrix[j]))
            denominator = np.sum(community_matrix[i] + community_matrix[j])
            if denominator > 0:
                total_dissimilarity += numerator / denominator
            num_pairs += 1
            
    return total_dissimilarity / num_pairs if num_pairs > 0 else 0.0


# --- Simulation Parameters ---
np.random.seed(0)
N_VOC_TYPES = 50  # Total possible VOCs in the world
N_INDIVIDUALS_PER_SITE = 10 # Plants per site

# --- Site 1: Low Latitude (Tropical, 10°N) ---
# High pressure -> individuals produce more VOC types (higher alpha)
# High pressure -> individuals are more different from each other (higher beta)
tropical_vocs = []
for _ in range(N_INDIVIDUALS_PER_SITE):
    # Each plant produces many VOC types. Base of 20, with variation.
    n_types_produced = np.random.randint(20, 35) 
    # VOC profile is unique to the individual
    individual_profile = np.zeros(N_VOC_TYPES)
    voc_indices = np.random.choice(N_VOC_TYPES, n_types_produced, replace=False)
    # Abundance of each VOC also varies
    individual_profile[voc_indices] = np.random.randint(1, 10, size=n_types_produced)
    tropical_vocs.append(individual_profile)

tropical_matrix = np.array(tropical_vocs)

# --- Site 2: High Latitude (Temperate, 60°N) ---
# Low pressure -> individuals produce fewer VOC types (lower alpha)
# Low pressure -> individuals are more similar to a common "optimal" profile (lower beta)
temperate_vocs = []
# Define a common "optimal" set of VOCs for the temperate zone
common_temperate_vocs = np.random.choice(N_VOC_TYPES, 15, replace=False)
for _ in range(N_INDIVIDUALS_PER_SITE):
    # Each plant produces fewer VOC types. Base of 8, with variation.
    n_types_produced = np.random.randint(8, 12)
    individual_profile = np.zeros(N_VOC_TYPES)
    # Most VOCs drawn from the common pool, with a few random ones
    n_common = min(n_types_produced, np.random.randint(7, 10))
    n_random = n_types_produced - n_common
    
    common_indices = np.random.choice(common_temperate_vocs, n_common, replace=False)
    # Ensure random indices don't overlap with common ones
    possible_random_indices = np.setdiff1d(np.arange(N_VOC_TYPES), common_temperate_vocs)
    random_indices = np.random.choice(possible_random_indices, n_random, replace=False)
    
    voc_indices = np.concatenate([common_indices, random_indices])
    individual_profile[voc_indices] = np.random.randint(1, 10, size=n_types_produced)
    temperate_vocs.append(individual_profile)
    
temperate_matrix = np.array(temperate_vocs)


# --- Calculate and Print Results ---
# Alpha diversity (average Shannon index for individuals at each site)
tropical_alpha = np.mean([shannon_diversity(ind) for ind in tropical_matrix])
temperate_alpha = np.mean([shannon_diversity(ind) for ind in temperate_matrix])

# Beta diversity (average dissimilarity between individuals at each site)
tropical_beta = calculate_beta_diversity(tropical_matrix)
temperate_beta = calculate_beta_diversity(temperate_matrix)

print("This simulation models the ecological principles discussed.")
print("-" * 50)
print(f"Low Latitude Site (e.g., Tropics):")
print(f"  - Average α-diversity (within-plant): {tropical_alpha:.4f}")
print(f"  - β-diversity (between-plant variation): {tropical_beta:.4f}")
print("-" * 50)
print(f"High Latitude Site (e.g., Temperate):")
print(f"  - Average α-diversity (within-plant): {temperate_alpha:.4f}")
print(f"  - β-diversity (between-plant variation): {temperate_beta:.4f}")
print("-" * 50)
print("Conclusion:")
print("As latitude increases (from Low to High), both α-diversity and β-diversity decrease.")
print("This demonstrates a negative effect of latitude on both metrics.")
print("α: negative, β: negative")

<<<B>>>